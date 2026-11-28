"""Select original outer blocks to refine at the same inner elastic width."""
import argparse
import gzip
import json
from math import fsum, isclose
from unittest.mock import patch

import early_pole_cold_screened_split_reference as split
from early_pole_prepared_hard_swap_trajectory_reference import run as prepared


def schedule(startup_outer_intervals,refined_blocks=None):
    if not 1<=startup_outer_intervals<=16:
        raise ValueError('Refine a proper count of original intervals.')
    if refined_blocks is None:
        return [1]*(4*startup_outer_intervals)+[4]*(16-startup_outer_intervals)
    if (not refined_blocks or len(set(refined_blocks))!=len(refined_blocks)
            or not all(type(i) is int and 1<=i<=16 for i in refined_blocks)):
        raise ValueError('Use distinct one-based original block indices from 1 to 16.')
    return [w for i in range(1,17) for w in ([1,1,1,1] if i in refined_blocks else [4])]


def evolve(bundle,case,startup_outer_intervals,refined_blocks=None):
    weights=schedule(startup_outer_intervals,refined_blocks)
    with gzip.open(bundle,'rt') as f:
        b=json.load(f)['branches']['steps']
    duration=b['duration']/b['steps']
    h=duration/64
    original=split.step
    runner=split.run
    def run(*a,**k):
        return runner(*a,**k,macro_step_weights=weights)
    def subcycle(model,y,dt,method='backward',energy_budget=None):
        q=round(dt/h)
        if q not in (1,4) or not isclose(dt,q*h,rel_tol=2e-15,abs_tol=0.):
            raise ValueError('Mismatched fixed inner width.')
        if energy_budget!=5e-15 or method!='backward':
            raise ValueError('Wrong individual elastic time profile.')
        out=[]
        for _ in range(q):
            y,d=original(model,y,dt/q,method,energy_budget=energy_budget)
            out.append(d)
        return y,dict(h=dt,elastic_substeps=out,
                      rho_drift=fsum(x['rho_drift'] for x in out),
                      entropy_gain=fsum(x['entropy_gain'] for x in out))
    with patch.object(split,'run',run),patch.object(split,'step',subcycle):
        x=prepared(bundle,case,len(weights),hold_screening=True)
    x['elastic_substep_energy_budget']=x.pop('elastic_energy_budget')
    if sum(len(r['elastic_substeps']) for r in x['elastic'])!=64:
        raise ValueError('Changed total elastic step count.')
    return x


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('bundle')
    p.add_argument('--case',required=True,choices=('native','local-daughter'))
    p.add_argument('--startup-outer-intervals',type=int,default=1)
    p.add_argument('--refined-blocks',type=int,nargs='+',
                   help='Override startup prefix: one-based original block indices.')
    a=p.parse_args()
    print(json.dumps(evolve(**vars(a)),indent=2,allow_nan=False))
