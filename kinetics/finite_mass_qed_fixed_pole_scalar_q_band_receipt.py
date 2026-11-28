"""Compare centered scalar-energy and number thresholds on complete q panels."""
import argparse
import gzip
import json
from math import fsum
from pathlib import Path

import numpy as np
from finite_mass_qed_timelike_pole_batched_absorption_reference import load,pole


def open_rows(path):
    if path.suffix=='.gz':
        with gzip.open(path,'rt') as s:return [json.loads(x) for x in s]
    return [json.loads(x) for x in path.read_text().splitlines()]


def collect(root,g,m,region,packed=False):
    w='.jsonl.gz' if packed else '.jsonl'
    v=open_rows(root/f'{g}-{region}-center-{m}{w}')
    if len(v)!=v[0]['total'] or [x['panel'] for x in v]!=list(range(len(v))):
        raise ArithmeticError(('Incomplete scalar q mesh',g,m,region))
    if any(a['limits'][1]!=b['limits'][0] for a,b in zip(v[:-1],v[1:],strict=True)):
        raise ArithmeticError('Scalar q panels overlap or leave holes')
    if any(x['grid']!=g or x['region']!=region for x in v):
        raise ArithmeticError('Unexpected scalar q scope')
    a=np.array([x['number_scalar_energy_and_centered_coefficients'] for x in v])
    if min(a[:,:,:2].ravel())<0.:
        raise ArithmeticError('Positive gain/loss sector turned negative')
    c=np.array([[fsum(a[:,j,k]) for k in range(3)] for j in range(3)])
    if np.max(abs((c[:,0]-c[:,1]-c[:,2])/np.maximum(c[:,0],c[:,1])))>3e-13:
        raise ArithmeticError('Separated statistical ledger is inconsistent')
    return c,[q['limits'] for q in v]


def case(root,g,region,packed=False):
    d,ch,_=load(g)
    M=d['checkpoints'][-1]['scalar_pole_over_T']
    pole_w=pole(ch,'T',.0015)[0]
    a,p=collect(root,g,0,region,packed)
    b,q=collect(root,g,1,region,packed)
    c,r=collect(root,g,2,region,packed)
    if p!=q or q!=r:
        raise ArithmeticError('Momentum refinement changed the exact topology')
    coarse=max(abs((a/b-1.).ravel()))
    mesh=max(abs((b/c-1.).ravel()))
    legacy=[]
    rel=[]
    for t in (a,b,c):
        if max(abs((M*t[0]+t[2])/t[1]-1.))>2e-14:
            raise ArithmeticError('Inconsistent full and centered scalar energy')
        z=np.longdouble(t[2])/np.longdouble(t[0])
        shift=(z[1]-z[2])/(M+z[2])
        n=t[0,1]/t[0,2]
        eta=t[1,1]/t[1,2]
        legacy.append([float(n),float(eta)])
        rel.append(float(shift))
        if abs(eta/n-1.-shift)>1e-14:
            raise ArithmeticError('Centered threshold drift cannot reproduce direct ratio')
    err=abs(rel[1]-rel[2])
    if mesh>1e-5 or (region=='cold' and not rel[2]>max(2e-9,40*err)):
        raise ArithmeticError(('Ambiguous band mesh or scalar heating sign',mesh,rel))
    return dict(grid=g,region=region,q_range=[p[0][0],p[-1][1]],panels=len(p),
                pole_k=.0015,pole_omega=pole_w,
                integrated_number_scalar_energy_and_excess_coefficients=c.tolist(),
                number_and_scalar_energy_threshold=legacy[2],
                relative_scalar_heat_neutral_shift=rel[2],mesh_discrepancy_of_shift=err,
                whole_coefficient_resolution_max_relative=float(mesh),
                coarse_to_intermediate_max_relative=float(coarse),all_centered_threshold_shifts=rel,
                common_prefactor_excluding_g2y2_and_dk=.0015**2/(64*np.pi**5))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    p.add_argument('--gzip-input',action='store_true')
    p.add_argument('--region',choices=('cold','shoulder'),default='cold')
    o=p.parse_args()
    print(json.dumps([case(o.root,g,o.region,o.gzip_input) for g in ('fine','third')],indent=2,allow_nan=False))
