"""Recover the initial-history-matched coarse pressure endpoint on the Gauss chart."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_startup_pilot import drive
from early_pole_endpoint_capacity_fast_step import step
from early_pole_endpoint_gauss_warm_chain_receipt import verify
from early_pole_endpoint_gauss_warm_comparison import gauss
from early_pole_endpoint_scalar_gauss2_live_fixture import ScalarGauss2
from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack
from early_pole_endpoint_stencil_remap_pilot import SPECS,load
from early_pole_resolved_parent_reference import S


def run(native,restart,ref_fine,old_report):
    lo,hi,n=.04,.05,1024
    orig,f=unpack(native,False,hi)
    mass=orig.state(hi,f)[0]
    d=ScalarGauss2(load('fine')[0],SPECS['fine'],2,4,1e-5,gauss_floor=.01)
    verify(old_report,ref_fine,d,lo,hi,n)
    if restart['x_end']!=lo or any(restart[k]!=ref_fine[k] for k in (
        'grid','charge_div','scalar_div','cut','gauss_floor','reference_x','nodes','scalar_cell_count','scalar_nodal_rule')):
        raise ArithmeticError('Inconsistent restart chart or clock')
    src=np.asarray(restart['nodal_occupation'])
    if src.shape!=(d.n,) or not np.all(np.isfinite(src)) or min(src)<0 or max(src[d.ns:])>1:
        raise ArithmeticError('Invalid shared restart state')
    trace,y=drive(d,hi,n,x_start=lo,initial=src,stepper=step)
    old=np.asarray([[z[k] for k in ('x','scalar_yield','continuity_residual','pair_number_error')]
                    for z in old_report['checkpoint_trajectories'][0]])
    now=np.asarray([[z[k] for k in ('x','scalar_yield','continuity_residual','pair_number_error')] for z in trace])
    gate=np.max(abs(old-now))
    if gate>5e-13 or np.max(abs(now[:,2]))>3e-10 or np.max(abs(now[:,3]))>3e-13:
        raise ArithmeticError('The recovered coarse branch does not reproduce its independent archive trace')
    u=fsum(d.ds.w*y[:d.ns])/S
    if abs(u/now[-1,1]-1)>4e-14:
        raise ArithmeticError('Recovered coarse budget failed')
    coarse=dict(ref_fine,nodal_occupation=y.tolist(),steps=n)
    a,b=[gauss(c,orig,mass,hi,expected_steps=k)[0] for c,k in ((coarse,n),(ref_fine,2*n))]
    return dict(x_start=lo,x_end=hi,steps=[n,2*n],fixed_pressure_mass=mass,maximum_trace_absolute_difference=gate,
                coarse_count_and_pressures=a.tolist(),fine_count_and_pressures=b.tolist(),
                coarse_over_fine_minus_one=(a/b-1).tolist(),
                scope='One conditional last-segment clock test with identical warm Gauss .04 restart; not a full-history error bound.'),coarse


if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('native',type=Path)
    a.add_argument('restart',type=Path)
    a.add_argument('fine_state',type=Path)
    a.add_argument('report',type=Path)
    a.add_argument('--coarse',type=Path,required=True)
    p=a.parse_args()
    r,c=run(*map(read,(p.native,p.restart,p.fine_state,p.report)))
    tmp=p.coarse.with_suffix(p.coarse.suffix+'.tmp')
    tmp.write_text(json.dumps(c,indent=2,allow_nan=False)+'\n')
    tmp.replace(p.coarse)
    print(json.dumps(r,indent=2,allow_nan=False))
