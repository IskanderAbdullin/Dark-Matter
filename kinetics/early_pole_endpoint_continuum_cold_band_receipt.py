"""Verify complete cold panels on an evolved endpoint snapshot."""
import argparse
import gzip
import json
from math import fsum
from pathlib import Path

import numpy as np


def read(path):
    op=gzip.open if path.suffix=='.gz' else open
    with op(path,'rt',encoding='utf-8') as f:
        return [json.loads(x) for x in f]


def reduce(path):
    x=read(path)
    if not x or any(r['total']!=len(x) for r in x) or [r['panel'] for r in x]!=list(range(len(x))):
        raise ArithmeticError('Incomplete cold interval')
    b=[r['bounds'] for r in x]
    if (b[0][0],b[-1][1])!=(1e-6,1.5e-5) or any(not lo<hi for lo,hi in b) or any(p[1]!=q[0] for p,q in zip(b[:-1],b[1:],strict=True)):
        raise ArithmeticError('Incomplete or overlapping retained support')
    tags=['grid','time_x','scalar_mass','omega','residue','mesh']
    if any(any(y[t]!=x[0][t] for t in tags) for y in x):
        raise ArithmeticError('Mixed endpoint or order')
    a=np.array([r['coeff'] for r in x])
    if a.shape!=(len(x),3,3) or not np.all(np.isfinite(a)) or min(x[0][t] for t in ('scalar_mass','omega','residue'))<=0:
        raise ArithmeticError('Invalid pole or coefficient card')
    if np.max(abs(a[:,:,0]-a[:,:,1]-a[:,:,2])/np.maximum(np.max(abs(a),axis=-1),1e-50))>4e-13:
        raise ArithmeticError('No corresponding affine ledger')
    z=np.array([[fsum(a[:,i,j]) for j in range(3)] for i in range(3)])
    if np.max(abs((x[0]['scalar_mass']*z[0]+z[2])/z[1]-1))>2e-14:
        raise ArithmeticError('Energy centering failed')
    return z,b,{t:x[0][t] for t in tags}


def run(files,archive):
    (a,ac,m),(b,bc,n)=[reduce(f) for f in files]
    if ac!=bc or any(m[k]!=n[k] for k in m if k!='mesh') or n['grid'] not in ('fine','third') or n['time_x']!=.05:
        raise ArithmeticError('Changed support or reference')
    if any(u>=v for u,v in zip(m['mesh'],n['mesh'],strict=True)) or min(np.min(a),np.min(b))<=0:
        raise ArithmeticError('Need positive moments on two distinct ordered meshes')
    err=np.max(abs(a/b-1))
    if err>6e-5:
        raise ArithmeticError(('Endpoint cold quadrature still underresolved',err))
    r={x['grid']:x for x in json.loads(archive.read_text())}[n['grid']]
    old=np.array(r['integrated_number_scalar_energy_and_excess_coefficients'])
    neutral=b[:2,1]/b[:2,2]
    split_a=(a[2,1]/a[0,1]-a[2,2]/a[0,2])/(n['scalar_mass']+a[2,2]/a[0,2])
    split_b=(b[2,1]/b[0,1]-b[2,2]/b[0,2])/(n['scalar_mass']+b[2,2]/b[0,2])
    if abs(split_b-split_a)>1e-10 or np.sign(split_a)!=np.sign(split_b):
        raise ArithmeticError('Paired neutrality split is unresolved')
    return dict(band=[1e-6,1.5e-5],full_pole_and_snapshot=n,accepted_meshes=[m['mesh'],n['mesh']],
                raw_number_scalar_energy_and_excess=b.tolist(),old_matched_pole_omega=r['pole_omega'],
                maximum_componentwise_mesh_difference=float(err),panels=len(bc),
                conditional_number_energy_neutral=neutral.tolist(),
                number_A_B_D_new_to_old=(b[0]/old[0]).tolist(),
                relative_neutrality_split=float(split_b),mesh_difference_in_relative_split=float(abs(split_b-split_a)),
                scope='One evolved snapshot and continuum closure, not spatial convergence or full reaction.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('raw_meshes',type=Path,nargs=2)
    p.add_argument('--archive',type=Path,default=Path(__file__).parent/'results'/'finite-mass-qed-fixed-pole-cold-scalar-band-2026-11-28.json')
    a=p.parse_args()
    print(json.dumps(run(a.raw_meshes,a.archive),indent=2,allow_nan=False))
