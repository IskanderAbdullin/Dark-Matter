"""Fixed T-pole q integration above the charged-transfer energy crossing."""
import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from finite_mass_qed_future_q_support_pilot import candidates
from finite_mass_qed_timelike_pole_angular_reference import angle
from finite_mass_qed_timelike_pole_batched_absorption_reference import load,pole


BANDS={'past_bridge':(.0003,.0015),'transfer_bridge':(.0015,.004),'first_future':(.004,.015),'opposite_pair':(.015,.1),'higher':(.1,1.)}
MESHES=((2,4,5,12),(4,8,8,18),(6,12,12,24))


def run(grid,mesh,band,out):
    d,ch,_=load(grid)
    m=d['checkpoints'][-1]['scalar_pole_over_T']
    k=.0015
    w,R=pole(ch,'T',k)
    lo,hi=BANDS[band]
    if band not in ('past_bridge','transfer_bridge') and not np.hypot(m,lo)>w:
        raise ArithmeticError('Unexpected transfer sign in future-band pilot')
    ss=np.asarray(d['scalar_momentum_edges'])
    occs=np.r_[0.,d['endpoint_scalar_occupations'],0.]
    edges,_=candidates(d,ch,k,w,m,lo,hi)
    n,cos,E,az=MESHES[mesh]
    xg,wg=roots_legendre(n)
    done=set()
    if out.exists():
        for line in out.read_text(encoding='utf-8').splitlines():
            obj=json.loads(line)
            j=obj['panel']
            if j in done or obj['grid']!=grid or obj['mesh']!=list(MESHES[mesh]) or obj['band']!=band or not (
                 0<=j<len(edges)-1 and obj['limits']==edges[j:j+2].tolist()):
                raise ArithmeticError(('Nonmatching partial receipt',j))
            done.add(j)
    with out.open('a',encoding='utf-8',buffering=1) as f:
        for j,(a,b) in enumerate(zip(edges[:-1],edges[1:],strict=True)):
            if j in done:
                continue
            nodes=(a+b)/2+(b-a)*xg/2
            wts=(b-a)*wg/2
            nmat=np.zeros((3,3))
            shell=[]
            for q,ww in zip(nodes,wts,strict=True):
                Es=np.hypot(m,q)
                ix=np.searchsorted(ss,q)
                shell.append(int(ix))
                occ=occs[ix]
                coeff,_,_,_=angle(ch,k,w,R,Es,'T',cos,(E,az),q,allow_future_spacelike=True)
                a0,b0,d0=coeff
                v0=ww*q*q/Es*np.array([(1.+occ)*a0,occ*b0,occ*d0+a0])
                nmat[0]+=v0
                nmat[1]+=Es*v0
                nmat[2]+=q*q/(Es+m)*v0
            if min(shell)!=max(shell):
                raise ArithmeticError(('Smearing exact archived shell',grid,band,j,shell))
            f.write(json.dumps(dict(grid=grid,mesh=MESHES[mesh],band=band,panel=j,total=len(edges)-1,
                                    limits=edges[j:j+2].tolist(),k=k,omega=w,
                                    number_scalar_energy_and_centered_coefficients=nmat.tolist()),
                               allow_nan=False)+'\n')
            f.flush()
            os.fsync(f.fileno())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=('fine','third'))
    p.add_argument('mesh',type=int,choices=range(len(MESHES)))
    p.add_argument('band',choices=BANDS)
    p.add_argument('out',type=Path)
    args=p.parse_args()
    run(args.grid,args.mesh,args.band,args.out)
