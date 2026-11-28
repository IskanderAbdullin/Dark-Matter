"""Recover the actual 1024-step mixed endpoint and compare with G7 on both clocks."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_startup_pilot import drive
from early_pole_endpoint_budgeted_warm_history import build
from early_pole_endpoint_capacity_fast_step import step
from early_pole_endpoint_gauss_warm_chain_receipt import verify_chart
from early_pole_endpoint_gauss_warm_comparison import gauss
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack
from early_pole_resolved_parent_reference import S


def run(native,prefix,control,restart,report,final,g7pair):
    b,bhash=read(native)
    archive,_=read(prefix)
    z=archive['contrasts'][-1]
    if (len(archive['mixed_trace'])!=5 or (z['x'],z['common_native_sha256'],z['steps'])!=(.05,bhash,[1024,2048])):
        raise ArithmeticError('Need the accepted complete paired histories and native measurement mass')
    d,meta=build(control)
    if meta['budget_reference_sha256']!=archive['budget_reference_sha256']:
        raise ArithmeticError('Original frozen budget selection changed')
    (s,shash),(r,rhash),(fin,fhash)=[read(p) for p in (restart,report,final)]
    prev,last=archive['mixed_trace'][-2:]
    if (shash,rhash,fhash)!=(prev['state_sha256'],last['report_sha256'],last['state_sha256']):
        raise ArithmeticError('Coarse replay cannot change the actual archived prefix')
    if any(s.get(k)!=v for k,v in meta.items()) or any(u['actual_previous_raw_sha256']!=shash for u in (r,fin)):
        raise ArithmeticError('Clock replay does not bind the exact actual layered chart or restart')
    old=np.asarray(s['nodal_occupation'])
    if (old.shape!=(d.n,) or not np.all(np.isfinite(old)) or min(old)<0 or max(old[d.ns:])>1
            or (s['x_end'],s['steps'])!=(.04,2048)):
        raise ArithmeticError('Invalid independently archived start')
    verify_chart(r,fin,d,.04,.05,1024,meta,d.repeats)
    b,w=unpack(b,False,.05)
    mass=b.state(.05,w)[0]
    if mass!=z['fixed_native_mass']:
        raise ArithmeticError('Did not reproduce the common measurement mass')
    trace,y=drive(d,.05,1024,x_start=.04,initial=old,stepper=step)
    words=('x','scalar_yield','continuity_residual','pair_number_error','max_charged','max_limited_extent')
    extracts=lambda p:np.asarray([[x[k] for k in words] for x in p])
    delta=np.max(abs(extracts(trace)-extracts(r['checkpoint_trajectories'][0])))
    if delta>5e-13 or abs((fsum(d.ds.w*y[:d.ns])/S)/trace[-1]['scalar_yield']-1)>4e-14:
        raise ArithmeticError('The recovered run does not match the independently archived original coarse run')
    co=dict(fin,steps=1024,nodal_occupation=y.tolist())
    a=np.array([gauss(c,b,mass,.05,expected_steps=k,budget_control=control)[0] for c,k in ((co,1024),(fin,2048))])
    g,gh=read(g7pair)
    if ((g['steps'],g['x'],g['measurement_mass'])!=([1024,2048],.05,mass)
            or g['charts'][1]['order']!=7 or g['charts'][1]['fine_state_sha256']!=archive['uniform_trace'][-1]['uncompressed_state_sha256']
            or g['charts'][1]['source_report_sha256']!=archive['uniform_trace'][-1]['source_report_sha256']):
        raise ArithmeticError('G7 clocks do not bind the complete original history')
    den=np.array(g['charts'][1]['values'])
    t=a/den-1
    if (np.max(abs(t[-1]-z['mixed255_over_G7_minus_one']))>5e-15 or
        np.max(abs(a[1]/z['moments'][0]-1))>5e-15):
        raise ArithmeticError('Did not reproduce independent fine contrast')
    return dict(x=.05,starts=.04,steps=[1024,2048],measurement_mass=mass,g7_clock_sha256=gh,
                mixed_report_sha256=rhash,mixed_fine_state_sha256=fhash,mixed_start_sha256=shash,
                mixed_coarse_trace_maximum_defect=delta,aligned_fields=['count','pressure','warm_pressure'],
                mixed_coarse_fine_moments=a.tolist(),mixed255_over_G7_minus_one=t.tolist(),
                fine_clock_minus_coarse=(t[1]-t[0]).tolist(),
                scope='Recovered same 1024-step branch from the same .04 state; no full-history error bound.'),co


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('native','prefix','control','restart','report','final','g7pair','coarse'):
        p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args()
    r,c=run(a.native,a.prefix,a.control,a.restart,a.report,a.final,a.g7pair)
    temp=a.coarse.with_suffix('.tmp')
    temp.write_text(json.dumps(c,indent=2,allow_nan=False)+'\n')
    temp.replace(a.coarse)
    print(json.dumps(r,indent=2,allow_nan=False))
