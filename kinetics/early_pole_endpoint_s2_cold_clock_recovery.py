"""Recover the missing paired coarse state for the S2 cold-band clock check."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_startup_pilot import drive
from early_pole_endpoint_capacity_fast_step import step
from early_pole_endpoint_charge_cold_clock import optical
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart
from early_pole_endpoint_stencil_remap_pilot import SPECS,X,load
from early_pole_resolved_parent_reference import S


def run(start,finish,report):
    a,sa=start; b,sb=finish; r,sr=report
    if (a['x_end'],b['x_end'],a['steps'],b['steps'],r['starting_x'],r['x_end'],r['coarse_steps'])!=(
            .04,.05,2048,2048,.04,.05,1024):
        raise ArithmeticError('Expected the exact stored final clock pair')
    d=ScalarSubshellChart(load('fine')[0],SPECS['fine'],2,2,1e-5)
    key=dict(grid='fine',charge_div=2,scalar_div=2,cut=1e-5,counts=[d.ns,d.nf,d.np],reference_x=X)
    if any(any(t[k]!=v for k,v in key.items()) or 'warm_floor' in t for t in (a,b,r)):
        raise ArithmeticError('The S2 clock chart changed')
    u,v=[np.asarray(t['endpoint_occupation']) for t in (a,b)]
    if any(z.shape!=(d.n,) or not np.all(np.isfinite(z)) or min(z)<0 or max(z[d.ns:])>1 for z in (u,v)):
        raise ArithmeticError('Invalid S2 occupations')
    trace,low=drive(d,.05,1024,x_start=.04,initial=u,stepper=step)
    cols=('x','scalar_yield','pair_number_error','continuity_residual','max_charged')
    t0,t1=r['checkpoint_trajectories']
    U,V=[np.asarray([[t[k] for k in cols] for t in h]) for h in (trace,t0)]
    high=fsum(d.ds.w*v[:d.ns])/S
    if U.shape!=(8,5) or V.shape!=U.shape or not np.all(np.isfinite(U)) or not np.all(np.isfinite(V)):
        raise ArithmeticError('Incomplete S2 clock trace')
    err=np.max(abs(U-V))
    if err>5e-13 or abs(high/t1[-1]['scalar_yield']-1)>4e-14:
        raise ArithmeticError('S2 recovery did not reconstruct the archived trajectory trace')
    h=[optical(d,z,.05) for z in (low,v)]
    return dict(x=.05,steps=[1024,2048],sha256_source=[sa,sb,sr],max_trace_defect=err,
                scalar_yield_coarse_fine=[trace[-1]['scalar_yield'],high],
                optical_parameters_coarse_fine=h,
                scope='Only a same-clock S2 endpoint recovery, not cold-band integration.'),dict(b,steps=1024,endpoint_occupation=low.tolist())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('start','fine','report','output'):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    r,s=run(*map(read,(a.start,a.fine,a.report)))
    temp=a.output.with_suffix('.tmp')
    temp.write_text(json.dumps(s,indent=2,allow_nan=False)+'\n')
    temp.replace(a.output)
    print(json.dumps(r,indent=2,allow_nan=False))
