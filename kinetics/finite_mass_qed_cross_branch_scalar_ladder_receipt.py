"""Verify full fixed-k scalar-momentum ladders while retaining mixed gate provenance."""
import argparse
import gzip
import json
from math import isfinite
from pathlib import Path

import numpy as np


MODES=(('T',.0015),('T',.00256),('L',.0012),('L',.00256))
Q=(0.,2.05e-6,4.1e-6,8.2e-6,4e-5,1e-4,3e-4,8e-4,.0012)
SUFFIX=('','-addendum','-selected-refinement')


def read(p):
    op=gzip.open if p.suffix=='.gz' else open
    with op(p,'rt',encoding='utf-8') as z:
        return [json.loads(s) for s in z]


def validate(row):
    a=np.asarray(row['raw_A_B_D'])
    if a.shape!=(3,) or not np.all(np.isfinite(a)) or min(a[:2])<0 or (
        abs(a[0]-a[1]-a[2])>5e-13*max(abs(a))):
        raise ArithmeticError('Invalid separated affine rates')
    if row['E_scalar']>=row['omega'] or not all(isfinite(row[k]) for k in ('E_scalar','omega')):
        raise ArithmeticError('These q samples require past transfer')
    if 'attempts' in row:
        arms=row['attempts']
        v,w=(np.asarray(d['raw']) for d in arms[-2:])
        diff=abs(v-w)
        err=diff/np.maximum(np.maximum(abs(v),abs(w)),1e-25)
        if not row['resolved_componentwise'] or w.tolist()!=a.tolist() or max(err)>=3e-6:
            raise ArithmeticError('Failed strict selected refinement')
        kind='componentwise'
    else:
        val=row['mesh_defect']
        if not 0<=val<=3e-6:
            raise ArithmeticError('Failed largest-component-scaled coarse comparison')
        diff=np.full(3,val*max(max(abs(a)),1e-25))
        kind='scaled_to_largest_component'
    f=row['scalar_occup']
    root=row['isotropic_conditional_population']
    expect=a[1]/(a[2]+a[0]/f) if f and a[2]+a[0]/f>0 and a[1]>0 else None
    if ((root is None)!=(expect is None) or (root is not None and abs(root/expect-1)>2e-14)):
        raise ArithmeticError('Conditional scalar source threshold mismatch')
    return diff,kind


def run(root,packed=True):
    pr='jsonl.gz' if packed else 'jsonl'
    v=[]
    for g in ('fine','third'):
        a={}
        for suf in SUFFIX:
            for r in read(root/f'q-ladder-{g}{suf}.{pr}'):
                key=(r['pol'],r['k'],r['q'])
                if r['grid']!=g or key in a:
                    raise ArithmeticError('Duplicate or wrong frozen archive entry')
                a[key]=r
        if set(a)!={(p,k,q) for p,k in MODES for q in Q}:
            raise ArithmeticError('The ladder is not a complete sampled set')
        drift={k:validate(r) for k,r in a.items()}
        for p,k in MODES:
            rows=[a[p,k,q] for q in Q]
            freq={r['omega'] for r in rows}
            if len(freq)!=1 or rows[0]['scalar_occup'] is not None or rows[0]['at_rest_angle_defect']>5e-13:
                raise ArithmeticError('Fixed mode momentum or scalar rest control differs')
            tests=[]
            for q,sign in ((4.1e-6,1),(8e-4,-1)):
                h=a[p,k,q]
                A,B,D=h['raw_A_B_D']
                f=h['scalar_occup']
                e=drift[p,k,q][0]
                n=1e-4
                signal=n*(D+A/f)-B
                separation=n*(e[2]+e[0]/f)+e[1]
                if signal*sign<=8*separation:
                    raise ArithmeticError('Observed mesh shift does not separate source sign at trial')
                tests.append(dict(q=q,threshold=h['isotropic_conditional_population'],scaled_net=float(signal),
                                  mesh_discrepancy_proxy=float(separation)))
            v.append(dict(grid=g,polarization=p,k=k,raw_q_samples=rows,
                          comparisons=[drift[p,k,q][1] for q in Q],common_trial_n=1e-4,two_mode_signs=tests))
    return dict(cases=v,nominal_common_opposite_sign_interval=[max(a['two_mode_signs'][0]['threshold'] for a in v),
                min(a['two_mode_signs'][1]['threshold'] for a in v)],
                scope='Discrete q probes at fixed poles. Observed mesh discrepancies are not rigorous error bars.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path,nargs='?',default=Path(__file__).parent/'results')
    p.add_argument('--raw',action='store_true')
    a=p.parse_args()
    print(json.dumps(run(a.root,not a.raw),indent=2,allow_nan=False))
