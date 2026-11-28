"""Carry a validated earlier coarse branch to the same fine final clock."""
import argparse
import hashlib
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


def rows(t):
    a=np.asarray([[r[k] for k in ('x','scalar_yield','pair_number_error','continuity_residual','max_charged','max_limited_extent')] for r in t])
    if (a.shape!=(8,6) or not np.all(np.isfinite(a)) or np.any(a[:,1:3]<0)
            or np.max(a[:,2])>3e-13 or np.max(abs(a[:,3]))>3e-10
            or np.any(a[:,4:]<0) or np.any(a[:,4:]>1)):
        raise ArithmeticError('Incomplete conservation, positivity or controller ledger')
    return a


def run(data,n):
    if len(data)!=5 or n not in (2,4):
        raise ValueError('Require independent S2/S4 records at .02, .04 and .05 and their paired reports')
    (a,ah),(b,bh),(c,ch),(r,rh),(s,sh)=data
    d=ScalarSubshellChart(load('fine')[0],SPECS['fine'],2,n,1e-5)
    key=dict(grid='fine',charge_div=2,scalar_div=n,cut=1e-5,counts=[d.ns,d.nf,d.np],reference_x=X)
    if any(any(z.get(k)!=v for k,v in key.items()) or z.get('warm_floor') is not None for z,_ in data):
        raise ArithmeticError('Early clock recovery would change geometry or chart')
    if any((p['x_end'],p['steps'])!=(x,2048) for p,x in zip((a,b,c),(.02,.04,.05),strict=True)):
        raise ArithmeticError('Unexpected archived fine interval')
    if any((p['starting_x'],p['x_end'],p['coarse_steps'])!=(lo,hi,1024)
           for p,lo,hi in [(r,.02,.04),(s,.04,.05)]):
        raise ArithmeticError('Archived clock reports do not bind the two prior windows')
    init,mid,ref=[np.asarray(p['endpoint_occupation']) for p in (a,b,c)]
    if any(z.shape!=(d.n,) or not np.all(np.isfinite(z)) or min(z)<0 or max(z[d.ns:])>1 for z in (init,mid,ref)):
        raise ArithmeticError('Invalid archived occupation')
    before=[rows(z) for z in r['checkpoint_trajectories']]
    after=[rows(z) for z in s['checkpoint_trajectories']]
    def count(y):
        return fsum(d.ds.w*y[:d.ns])/S
    for y,t in ((mid,before[1]),(ref,after[1])):
        if abs(count(y)/t[-1,1]-1)>4e-14 or abs(max(y[d.ns:])-t[-1,4])>3e-14:
            raise ArithmeticError('Archived reference states disagree with fine traces')
    u,A=drive(d,.04,1024,x_start=.02,initial=init,stepper=step)
    err=float(np.max(abs(rows(u)-before[0])))
    if err>5e-13 or abs(count(A)/u[-1]['scalar_yield']-1)>4e-14:
        raise ArithmeticError('Recovered prior coarse state does not match original prior coarse trace')
    old=dict(b,steps=1024,endpoint_occupation=A.tolist(),prior_clock_control='coarse_02_to_04')
    wire=json.dumps(old,indent=2,allow_nan=False)+'\n'
    p,B=drive(d,.05,2048,x_start=.04,initial=A,stepper=step)
    if rows(p)[-1,1]<=0 or abs(count(B)/p[-1]['scalar_yield']-1)>4e-14 or min(B)<0 or max(B[d.ns:])>1:
        raise ArithmeticError('Alternative terminal ledger does not close')
    final=dict(c,endpoint_occupation=B.tolist(),prior_clock_control='coarse_02_to_04',
               actual_prior_state_sha256=hashlib.sha256(wire.encode()).hexdigest())
    poles=[optical(d,y,.05) for y in (B,ref)]
    return dict(x=.05,grid='fine',scalar_div=n,charge_div=2,inputs_sha256=[ah,bh,ch,rh,sh],
                alternate_path=[1024,2048],original_path=[2048,2048],archive_prior_coarse_trace_defect=err,
                earlier_coarse_recovery=u,alternative_final_trace=p,
                final_alternate_and_original_yield=[count(B),count(ref)],
                final_alternate_and_original_poles=poles,
                occupation_l1_difference_relative_to_reference=fsum(d.ds.w*abs(B[:d.ns]-ref[:d.ns]))/(S*count(ref)),
                scope='Coarsens only the .02→.04 clock; same fine final clock. Not a full-history extrapolation or physical prediction.'),wire,final


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs',type=Path,nargs=5,required=True,metavar=('X02','X04','X05','X04_REPORT','X05_REPORT'))
    p.add_argument('--level',type=int,choices=(2,4),required=True)
    p.add_argument('--earlier',type=Path,required=True)
    p.add_argument('--terminal',type=Path,required=True)
    a=p.parse_args()
    r,u,v=run([read(p) for p in a.inputs],a.level)
    for name,t in ((a.earlier,u),(a.terminal,json.dumps(v,indent=2,allow_nan=False)+'\n')):
        h=name.with_suffix('.part')
        h.write_text(t)
        h.replace(name)
    print(json.dumps(r,indent=2,allow_nan=False))
