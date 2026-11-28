"""Pair held-mass warm Gauss2 and Gauss3 numerical histories and frozen limit moments."""
import argparse
import json
from pathlib import Path

import numpy as np

from early_pole_endpoint_scalar_subshell_spectrum_audit import read


ORDER=('count','pressure','pressure_q_ge_.01')


def match(d2,d3):
    if len(d2)!=4 or len(d3)!=4:
        raise ArithmeticError('Need four complete matched comparison states')
    out=[]
    for x,a,b in zip((.01,.02,.04,.05),d2,d3,strict=True):
        if (a['x'],b['x'])!=(x,x) or a['fixed_mass']!=b['fixed_mass'] or (
                tuple(a['columns'])!=ORDER or tuple(b['columns'])!=ORDER
                or a['S4_reference']!=b['S4_reference'] or a.get('warm_gauss_order',2)!=2 or b.get('warm_gauss_order')!=3):
            raise ArithmeticError('Did not compare the same reference states, quantities and mass')
        s,t=[np.array(c['warm_gauss_moments']) for c in (a,b)]
        lim=[np.array(c['warm_gauss_conditional_pressure_envelope']['pressure_total_min_max']) for c in (a,b)]
        if min(s)<=0 or min(t)<=0 or any(not 0<v[0]<=v[1] for v in lim):
            raise ArithmeticError('Nonpositive conditional moments')
        if any(not lo<=q<=hi for (lo,hi),q in zip(lim,(s[1],t[1]),strict=True)):
            raise ArithmeticError('Positive representative left its moment interval')
        j=[lim[1][0]/lim[0][1]-1,lim[1][1]/lim[0][0]-1]
        u,v=[{r['warm_div']:r['moments'] for r in c['geometric_warm_charts']} for c in (a,b)]
        if u!=v or set(u)!={2,4}:
            raise ArithmeticError('Both comparisons must include the actual common W2 and W4')
        out.append(dict(x=x,measurement_mass=a['fixed_mass'],Gauss3_over_Gauss2_minus_one=(t/s-1).tolist(),
                        pressure_contrast_interval_at_discrete_cell_count_and_energy=j,
                        Gauss3_over_W4_minus_one=(t/np.array(u[4])-1).tolist()))
    return out


def stress(a,b):
    keys=('fixed_evolved_x','tail_cut','reference_fitted_mass','columns')
    if any(a[k]!=b[k] for k in keys) or a['frozen_limit_columns']!=b['frozen_limit_columns']:
        raise ArithmeticError('Mismatched frozen observation conventions')
    aa,bb=a['measurements'],b['measurements']
    if len(aa)!=len(bb):
        raise ArithmeticError('Mismatched frozen weight grids')
    out=[]
    for x,y in zip(aa,bb,strict=True):
        if x['test_mass']!=y['test_mass'] or x['rows'][:-1]!=y['rows'][:-1]:
            raise ArithmeticError('Frozen weights or shared geometric histories changed')
        u,v=x['rows'][-1],y['rows'][-1]
        if (u['warm_chart'],v['warm_chart'])!=('Gauss2','Gauss3'):
            raise ArithmeticError('Incorrect matched warm nodal checks')
        out.append(dict(test_mass=x['test_mass'],
                        Gauss3_over_Gauss2_minus_one=(np.array(v['totals'])/np.array(u['totals'])-1).tolist()))
    u,v=[z['limits'][-1]['finite_moment_limits'] for z in (a,b)]
    return dict(frozen_columns=a['columns'],rows=out,limit_columns=a['frozen_limit_columns'],
                Gauss3_over_Gauss2_limit_minus_one=(np.array(v)/np.array(u)-1).tolist())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--order2',type=Path,nargs=4,required=True)
    p.add_argument('--order3',type=Path,nargs=4,required=True)
    p.add_argument('--stress',type=Path,nargs=2,required=True)
    a=p.parse_args()
    print(json.dumps(dict(dynamic_rows=match([read(p) for p in a.order2],[read(p) for p in a.order3]),
                          frozen_stress=stress(*map(read,a.stress)),
                          scope='Two fully evolved numerical warm charts and their conditioned positive reconstructions; no momentum convergence or physical late-epoch prediction.'),indent=2,allow_nan=False))
