"""Resolve scalar warm-cell nodal biases on one unchanged frozen histogram."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_scalar_spectral_moment_budget import LABELS,integrate,load_card
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart
from early_pole_endpoint_stencil_remap_pilot import SPECS,load


def run(p):
    src,state,m=load_card(p,4)
    a=load('fine')[0]
    ref=integrate(src,state,m,48).sum(axis=0)
    out=[]
    former=None
    exact_support=src.ds.ed[src.ds.ed<=.01]
    for lev in (1,2,4,8):
        d=ScalarSubshellChart(a,SPECS['fine'],2,4,1e-5,warm_floor=.01,warm_div=lev)
        if any(not np.array_equal(getattr(src,k),getattr(d,k)) for k in ('qd','qp','wd','wp')) or (
                not np.array_equal(d.ds.ed[d.ds.ed<=.01],exact_support)):
            raise ArithmeticError('Warm-only residual changed cold or charge geometry')
        q=d.ds.p
        parent=np.searchsorted(src.ds.ed,q)-1
        f=np.r_[state[parent],state[src.ns:]]
        cells=integrate(d,f,m,48)
        if np.max(abs(cells.sum(axis=0)/ref-1))>6e-14:
            raise ArithmeticError('The original declared histogram changed after splitting')
        N=d.ds.w*f[:d.ns]
        E=np.hypot(q,m)
        nodal=[fsum(N*g) for g in (q**0,q,q*q,q*q/(E+m),q*q/E)]
        err=nodal/ref-1
        if former is not None and (abs(err[2])>=abs(former[2]) or abs(err[4])>=abs(former[4])):
            raise ArithmeticError('Warm pressure or second moment failed to improve')
        out.append(dict(warm_div=lev,scalar_bins=d.ns,relative_error_of_nodal_moment=err.tolist(),
                       invariance_of_exact_total=float(np.max(abs(cells.sum(axis=0)/ref-1)))))
        former=err
    return dict(cold_scalar_div=4,charge_div=2,warm_floor_for_complete_native_cells=.01,
                moment_columns=LABELS,cards=out,
                scope='One frozen piecewise-constant histogram; a representational audit, not re-evolution or dynamic convergence.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('state',type=Path)
    print(json.dumps(run(p.parse_args().state),indent=2,allow_nan=False))
