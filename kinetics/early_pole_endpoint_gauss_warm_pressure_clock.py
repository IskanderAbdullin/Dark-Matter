"""Recover the initial-history-matched coarse pressure endpoint on the Gauss chart."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_startup_pilot import drive
from early_pole_endpoint_capacity_fast_step import step
from early_pole_endpoint_gauss_warm_chain_receipt import verify,read as read_and_hash
from early_pole_endpoint_gauss_warm_comparison import chart_for,gauss
from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack
from early_pole_resolved_parent_reference import S


def run(native,restart,ref_fine,old_report,*,first_continuation=False,restart_sha=None):
    lo,hi,n=(.0075,.01,512) if first_continuation else (.04,.05,1024)
    if first_continuation and (not restart_sha or restart.get('steps')!=4096
                              or any(q.get('actual_previous_raw_sha256')!=restart_sha for q in (ref_fine,old_report))):
        raise ArithmeticError('First continuation does not bind the original startup checkpoint')
    orig,f=unpack(native,False,hi)
    mass=orig.state(hi,f)[0]
    order,d,_=chart_for(ref_fine)
    verify(old_report,ref_fine,d,lo,hi,n,order=None if order==2 else order)
    if order!=2 and restart.get('warm_gauss_order')!=order:
        raise ArithmeticError('Restarts cannot change the warm Gauss order')
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
    kind={} if order==2 else dict(warm_gauss_order=order)
    return dict(x_start=lo,x_end=hi,steps=[n,2*n],fixed_pressure_mass=mass,maximum_trace_absolute_difference=gate,**kind,
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
    a.add_argument('--first-continuation',action='store_true',help='replay .0075 to .01 instead of .04 to .05')
    p=a.parse_args()
    k,l,m,n=map(read,(p.native,p.restart,p.fine_state,p.report))
    digest=read_and_hash(p.restart)[1] if p.first_continuation else None
    r,c=run(k,l,m,n,first_continuation=p.first_continuation,restart_sha=digest)
    tmp=p.coarse.with_suffix(p.coarse.suffix+'.tmp')
    tmp.write_text(json.dumps(c,indent=2,allow_nan=False)+'\n')
    tmp.replace(p.coarse)
    print(json.dumps(r,indent=2,allow_nan=False))
