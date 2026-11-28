"""Validate the missing scalar interval before joining covered fixed-pole bands."""
import argparse
import gzip
import json
from math import fsum
from pathlib import Path

import numpy as np


COLD=(1e-6,1.5e-5)
GAP=(1.5e-5,4e-5)
SHOULDER=(4e-5,3e-4)


def read(path):
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'rt',encoding='utf-8') as f:
        return [json.loads(line) for line in f]


def aggregate(rows,grid,band,bounds,key):
    indices=[r['panel'] for r in rows]
    if indices!=list(range(len(rows))) or not len(rows) or any(
        r['total']!=len(rows) or r['grid']!=grid or r['region']!=band or r['mesh']!=rows[0]['mesh']
        for r in rows):
        raise ArithmeticError(('Incomplete or mixed mesh',grid,band,len(rows)))
    supports=[r['limits'] for r in rows]
    if (supports[0][0],supports[-1][1])!=bounds or any(
        x[1]!=y[0] for x,y in zip(supports[:-1],supports[1:],strict=True)):
        raise ArithmeticError(('Noncontiguous coverage',grid,band))
    A=np.array([r[key] for r in rows])
    if not np.all(np.isfinite(A)) or np.min(A[:,:,:2])<0:
        raise ArithmeticError(('Non-positive input',grid,band))
    err=abs(A[:,:,0]-A[:,:,1]-A[:,:,2])/np.maximum(np.max(abs(A),axis=2),np.finfo(float).tiny)
    if np.max(err)>3e-13:
        raise ArithmeticError(('Forward-inverse closure',grid,band,np.max(err)))
    total=np.array([[fsum(A[:,i,j]) for j in range(3)] for i in range(A.shape[1])])
    return total,supports


def call(files,grid,band,span,key):
    a,b=[aggregate(read(f),grid,band,span,key) for f in files]
    if a[1]!=b[1]:
        raise ArithmeticError(('Disagreement on a shell edge',grid,band))
    r,high=a[0],b[0]
    defect=np.max(abs(r-high)/np.maximum(abs(high),np.finfo(float).tiny))
    if defect>6e-5:
        raise ArithmeticError(('Unresolved whole coefficient',grid,band,defect))
    T=high[:2,1]/high[:2,2]
    if not np.all(np.isfinite(T)) or np.min(T)<=0:
        raise ArithmeticError(('Unsupported conditional neutrality',grid,band,T))
    neutral_L=r[:2,1]/r[:2,2]
    return dict(panel_count=len(b[1]),support=span,raw=high.tolist(),
                neutral_n=T.tolist(),neutral_mesh_drift=abs(neutral_L-T).tolist(),
                max_coefficient_mesh_defect=float(defect))


def run(shared,hole,pack):
    if COLD[1]!=GAP[0] or GAP[1]!=SHOULDER[0]:
        raise ArithmeticError('Gap is not the exact complement of the named intervals')
    output=[]
    for g in ('fine','third'):
        configs=((shared,'cold',COLD,'number_scalar_energy_and_centered_coefficients','-center-',(1,2),'.gz'),
                 (hole,'gap',GAP,'number_scalar_energy_and_centered_coefficients','-center-',(1,2),pack),
                 (shared,'shoulder',SHOULDER,'number_and_scalar_energy_coefficients','-', (0,1),'.gz'))
        pieces=[]
        for d,name,span,key,pattern,meshes,suffix in configs:
            t=call([d/f'{g}-{name}{pattern}{m}.jsonl{suffix}' for m in meshes],g,name,span,key)
            pieces.append(dict(name=name,**t))
        coef=np.array([x['raw'][:2] for x in pieces])
        total=np.array([[fsum(coef[:,i,j]) for j in range(3)] for i in range(2)])
        neutral=total[:,1]/total[:,2]
        terms=coef[:,:,2]*neutral[None,:]-coef[:,:,1]
        skip=coef[0]+coef[2]
        no_gap=skip[:,1]/skip[:,2]
        relative_gap_shift=neutral/no_gap-1
        if np.max(abs(relative_gap_shift))>2e-6:
            raise ArithmeticError(('Gap is no longer negligible for these frozen archives',g,relative_gap_shift))
        if np.any(abs(np.sum(terms,axis=0))>2e-12*np.max(abs(terms),axis=0)):
            raise ArithmeticError(('Retained union does not close',g,terms))
        output.append(dict(grid=g,pieces=pieces,joined_support=(COLD[0],SHOULDER[1]),raw=total.tolist(),
                           neutral_number_and_energy=neutral.tolist(),
                           neutral_number_and_energy_skipping_gap=no_gap.tolist(),
                           relative_neutral_shift_including_gap=relative_gap_shift.tolist(),
                           scalar_energy_at_joined_number_neutral=float(neutral[0]*total[1,2]-total[1,1]),
                           corresponding_signed_transfer_by_piece=terms.tolist(),
                           at_n_1e_8=(coef[:,:,2]*1e-8-coef[:,:,1]).tolist()))
    return output


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    dest=Path(__file__).parent/'results'
    p.add_argument('--shared',type=Path,default=dest)
    p.add_argument('--gap',type=Path,default=dest)
    p.add_argument('--raw-gap',action='store_true')
    o=p.parse_args()
    print(json.dumps(run(o.shared,o.gap,'' if o.raw_gap else '.gz'),indent=2,allow_nan=False))
