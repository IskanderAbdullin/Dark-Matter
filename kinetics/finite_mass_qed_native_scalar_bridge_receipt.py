"""Require complete native bridge panels and compare all three energy ledgers."""
import argparse
import gzip
import json
from math import fsum
from pathlib import Path

import numpy as np


BANDS={'past-bridge':(.0003,.0015),'transfer-bridge':(.0015,.004)}
MESHES=((2,4,5,12),(4,8,8,18),(6,12,12,24))


def read(path):
    op=gzip.open if path.suffix=='.gz' else open
    with op(path,'rt',encoding='utf-8') as f:
        return [json.loads(s) for s in f]


def panels(path,grid,band,mesh,omega,m):
    y=read(path)
    if not y or [r['panel'] for r in y]!=list(range(len(y))):
        raise ArithmeticError('Missing, duplicate, or misordered raw panel')
    tags=dict(grid=grid,band=band.replace('-','_'),mesh=list(MESHES[mesh]),
              k=.0015,omega=omega,total=len(y))
    if any(any(q[key]!=v for key,v in tags.items()) for q in y):
        raise ArithmeticError('Wrong archive, pole, order, or incomplete raw file')
    lim=[r['limits'] for r in y]
    if (lim[0][0],lim[-1][1])!=BANDS[band] or any(not l<h for l,h in lim) or any(
            a[1]!=b[0] for a,b in zip(lim[:-1],lim[1:],strict=True)):
        raise ArithmeticError('Bridge has a gap, reversed panel, or changed limits')
    v=np.asarray([r['number_scalar_energy_and_centered_coefficients'] for r in y])
    if v.shape!=(len(y),3,3) or not np.all(np.isfinite(v)) or np.min(v[:,:,:2])<0:
        raise ArithmeticError('Invalid gross process coefficient')
    s=np.maximum(np.max(abs(v),axis=-1),1e-65)
    if np.max(abs(v[:,:,0]-v[:,:,1]-v[:,:,2])/s)>5e-13:
        raise ArithmeticError('Paired forward and inverse columns did not close')
    if np.max(abs(m*v[:,0]+v[:,2]-v[:,1])/s[:,1,None])>4e-13:
        raise ArithmeticError('Centered scalar heat did not close')
    t=np.hypot(lim,m)-omega
    if np.any(t[:,0]*t[:,1]<-1e-25):
        raise ArithmeticError('Transfer-energy crossing inside an uncut panel')
    parts=np.array([[fsum(v[:,i,j]) for j in range(3)] for i in range(3)])
    h=(omega-m)*v[:,0]-v[:,2]
    check=omega*v[:,0]-v[:,1]
    err=np.max(abs(h-check)/np.maximum(omega*s[:,0,None]+s[:,1,None],1e-65))
    if err>3e-13:
        raise ArithmeticError('Independent charged-energy ledger did not close')
    return lim,parts,np.array([fsum(h[:,j]) for j in range(3)]),t


def run(root,cold,band,mesh,pack,grids=('fine','third')):
    ref={u['grid']:u for u in json.loads(cold.read_text())}
    out=[]
    if tuple(grids) not in (('fine',),('third',),('fine','third')):
        raise ValueError('Select one valid grid or both in canonical order')
    for g in grids:
        r=ref[g]
        old=np.asarray(r['integrated_number_scalar_energy_and_excess_coefficients'])
        w=r['pole_omega']
        m=(old[1,0]-old[2,0])/old[0,0]
        a,b=[panels(root/f'{g}-{band}-{i}.jsonl{pack}',g,band,i,w,m) for i in mesh]
        if a[0]!=b[0]:
            raise ArithmeticError('Refinement has altered topology')
        test=lambda u,v: abs(u-v)/np.maximum(np.maximum(abs(u),abs(v)),1e-65)
        s= test(a[1],b[1])
        q= test(a[2],b[2])
        if max(np.max(s),np.max(q))>3e-5:
            raise ArithmeticError(('Unresolved bridge',g,float(np.max(s)),float(np.max(q))))
        # A root exists only if inverse and slope have the same sign.
        D=[b[1][0,2],b[1][1,2],b[2][2]]
        B=[b[1][0,1],b[1][1,1],b[2][1]]
        roots=[float(x/y) if x*y>0 else None for x,y in zip(B,D,strict=True)]
        out.append(dict(grid=g,band=band,scope=list(BANDS[band]),omega=w,m=m,accepted_orders=mesh,
                        panels=len(b[0]),past_or_future_counts=[int(np.sum(np.mean(b[3],axis=1)<0)),
                                                               int(np.sum(np.mean(b[3],axis=1)>0))],
                        number_scalar_energy_and_centered_coefficients=b[1].tolist(),
                        charged_energy_coefficients=b[2].tolist(),mesh_relative_drift=s.tolist(),
                        charged_mesh_relative_drift=q.tolist(),number_scalar_and_charged_neutral_n=roots))
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    d=Path(__file__).parent/'results'
    p.add_argument('band',choices=BANDS)
    p.add_argument('--root',type=Path,default=d)
    p.add_argument('--cold',type=Path,default=d/'finite-mass-qed-fixed-pole-cold-scalar-band-2026-11-28.json')
    p.add_argument('--meshes',nargs=2,type=int,choices=(0,1,2),default=(0,1))
    p.add_argument('--raw',action='store_true')
    p.add_argument('--grids',nargs='+',choices=('fine','third'),default=('fine','third'))
    a=p.parse_args()
    if a.meshes[0]>=a.meshes[1]:
        raise ValueError('Quadrature refinement must increase')
    print(json.dumps(run(a.root,a.cold,a.band,a.meshes,'' if a.raw else '.gz',a.grids),indent=2,allow_nan=False))
