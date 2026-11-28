"""Compare Gauss5/Gauss6 on two recovered first-continuation clocks."""
import argparse
import json
from pathlib import Path

import numpy as np


def values(replay,compare,n):
    if ((replay['x_start'],replay['x_end'],replay['steps'],replay['warm_gauss_order'])!=(.0075,.01,[512,1024],n)
        or (compare['x'],compare['warm_gauss_order'])!=(.01,n)
        or replay['fixed_pressure_mass']!=compare['fixed_mass']
        or not 0.<=replay['maximum_trace_absolute_difference']<=5e-13):
        raise ArithmeticError('Wrong paired first-restart epoch, warm order or fixed observable mass')
    v=np.array([replay['coarse_count_and_pressures'],replay['fine_count_and_pressures']])
    if v.shape!=(2,3) or not np.all(np.isfinite(v)) or v.min()<=0.:
        raise ArithmeticError('Invalid paired count and pressure observations')
    if (np.max(abs(v[0]/v[1]-1.-replay['coarse_over_fine_minus_one']))>1e-15
        or np.max(abs(v[1]/compare['warm_gauss_moments']-1.))>5e-15):
        raise ArithmeticError('Pressure replay disagrees with independent fine endpoint')
    return v


def run(replay,compare):
    a,b=[values(r,c,k) for r,c,k in zip(replay,compare,(5,6),strict=True)]
    if compare[0]['fixed_mass']!=compare[1]['fixed_mass'] or compare[0]['S4_reference']!=compare[1]['S4_reference']:
        raise ArithmeticError('First-restart charts did not use the same fixed S4 source')
    v=b/a-1.
    return dict(grid='fine',x_start=.0075,x_end=.01,steps=[512,1024],levels=[5,6],
                fixed_mass=compare[0]['fixed_mass'],fields=['count','pressure','warm_pressure'],
                G6_over_G5_minus_one_at_each_clock=v.tolist(),
                fine_minus_coarse_of_G6_over_G5=(v[1]-v[0]).tolist(),
                G5_and_G6_individual_coarse_over_fine_minus_one=[(m[0]/m[1]-1.).tolist() for m in (a,b)],
                scope='Each order holds its own earlier startup state fixed. This does not test startup pressure error or subsequent evolution.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('replay','comparison'):
        p.add_argument('--'+k,type=Path,nargs=2,required=True)
    a=p.parse_args()
    def cards(paths):
        return [json.loads(x.read_text()) for x in paths]
    print(json.dumps(run(cards(a.replay),cards(a.comparison)),indent=2,allow_nan=False))
