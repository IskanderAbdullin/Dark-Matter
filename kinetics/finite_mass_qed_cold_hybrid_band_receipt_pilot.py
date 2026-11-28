"""Reduce complete hybrid cold bands and attribute their positive integrals."""
import argparse
import gzip
import json
from itertools import permutations
from math import fsum
from pathlib import Path

import numpy as np


GRIDS=('fine','third')


def collect(path,ch,geo):
    f=gzip.open if path.suffix=='.gz' else open
    with f(path,'rt',encoding='utf-8') as z:
        data=[json.loads(s) for s in z]
    if not data or len(data)!=data[0]['total'] or [r['panel'] for r in data]!=list(range(len(data))):
        raise ArithmeticError(('Missing hybrid panels',path,len(data)))
    if any(r['charged_grid']!=ch or r['dispersion_grid']!=geo or r['scalar_grids']!=list(GRIDS) or
           any(r[k]!=data[0][k] for k in ('total','mesh','k','w','m_scalar')) for r in data):
        raise ArithmeticError('Wrong archive/dispersion/mesh')
    bounds=[r['bounds'] for r in data]
    if (bounds[0][0],bounds[-1][1])!=(1e-6,1.5e-5) or any(
        a[1]!=b[0] for a,b in zip(bounds[:-1],bounds[1:],strict=True)):
        raise ArithmeticError('Cold support not fully covered')
    a=np.array([r['number_scalar_energy_excess'] for r in data])
    if a.shape!=(len(data),2,3,3) or np.min(a[:,:,:,:2])<0 or not np.all(np.isfinite(a)):
        raise ArithmeticError('Invalid gross integral')
    if np.max(abs(a[:,:,:,0]-a[:,:,:,1]-a[:,:,:,2])/np.maximum(np.max(abs(a),axis=3),1e-55))>4e-13:
        raise ArithmeticError('Separate inverse and forward do not close')
    result=np.array([[[fsum(a[:,j,k,i]) for i in range(3)] for k in range(3)] for j in range(2)])
    reconstructed=data[0]['m_scalar']*result[:,0]+result[:,2]
    if np.max(abs((reconstructed-result[:,1])/result[:,1]))>3e-14:
        raise ArithmeticError('Centered/full energy moments do not agree')
    return result,bounds


def effects(v):
    e=np.zeros(3)
    for order in permutations(range(3)):
        g=['third']*3
        a=np.log(v[tuple(g)])
        for i in order:
            g[i]='fine'
            b=np.log(v[tuple(g)])
            e[i]+=(b-a)/6
            a=b
    audit=np.log(v[('fine',)*3]/v[('third',)*3])
    if abs(sum(e)-audit)>2e-13:
        raise ArithmeticError('Ordering-average decomposition does not close')
    return e.tolist()


def run(root,native,meshes,pack):
    reference={r['grid']:np.asarray(r['integrated_number_scalar_energy_and_excess_coefficients'])
               for r in json.loads(native.read_text())}
    contributions={}
    rows=[]
    for ch in GRIDS:
        for geo in GRIDS:
            a,b=[collect(root/f'{ch}-{geo}-{mesh}.jsonl{pack}',ch,geo) for mesh in meshes]
            if a[1]!=b[1]:
                raise ArithmeticError('Hybrid support changed under mesh refinement')
            drift=float(np.max(abs((a[0]-b[0])/b[0])))
            if drift>6e-5 or np.min(b[0])<=0:
                raise ArithmeticError(('Unresolved coefficients',ch,geo,drift))
            for i,g in enumerate(GRIDS):
                p=b[0][i]
                nt=p[:2,1]/p[:2,2]
                entry=dict(charged_grid=ch,dispersion_grid=geo,scalar_histogram=g,
                           raw_number_scalar_energy_excess=p.tolist(),conditional_neutral_n=nt.tolist(),
                           mesh_max_discrepancy=drift,q_panels=len(b[1]))
                if ch==geo==g:
                    e=float(np.max(abs((p-reference[g])/reference[g])))
                    if e>1e-5:
                        raise ArithmeticError(('Reference main path mismatch',g,e))
                    entry['native_relative_replay_defect']=e
                contributions[ch,geo,g]=(p[0,0],p[0,1],p[0,2],nt[0])
                rows.append(entry)
    return dict(scope='Controlled substitutions of archived cells within one cold band and fixed photon k.',
                factor_order=['charged_histogram','dispersion_pair','scalar_histogram'],
                quantity_order=['forward_coeff','inverse_coeff','affine_coeff','conditional_number_neutral'],
                log_effects=[effects({p:x[i] for p,x in contributions.items()}) for i in range(4)],
                data=rows)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',type=Path)
    parser.add_argument('native',type=Path)
    parser.add_argument('--meshes',type=int,nargs=2,default=[0,1])
    parser.add_argument('--raw',action='store_true')
    inp=parser.parse_args()
    print(json.dumps(run(inp.root,inp.native,inp.meshes,'' if inp.raw else '.gz'),indent=2,allow_nan=False))
