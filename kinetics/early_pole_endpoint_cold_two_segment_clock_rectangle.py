"""Resolve the missing combined coarse corner and check two cold momentum points."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_startup_pilot import drive
from early_pole_endpoint_capacity_fast_step import step
from early_pole_endpoint_charge_cold_clock import optical
from early_pole_endpoint_charge_cold_clock_points import setup,point
from early_pole_endpoint_cold_prior_clock_restart import rows
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart
from early_pole_endpoint_stencil_remap_pilot import SPECS,X,load
from early_pole_resolved_parent_reference import S


def decompose(value):
    v=np.asarray(value)
    if v.shape[:2]!=(2,2) or np.min(v)<=0 or not np.all(np.isfinite(v)):
        raise ArithmeticError('The two-clock grid requires positive variables on all four corners')
    h=np.log(v)
    p=h[1,0]-h[0,0]
    f=h[0,1]-h[0,0]
    mix=h[1,1]-h[1,0]-h[0,1]+h[0,0]
    if np.max(abs(p+f+mix-h[1,1]+h[0,0]))>3e-14:
        raise ArithmeticError('Clock interaction decomposition does not reconstruct the observed combined shift')
    return dict(corner_values=v.tolist(),earlier_only_log_shift=p.tolist(),final_only_log_shift=f.tolist(),
                interaction_log_shift=mix.tolist(),both_minus_original_log_shift=(h[1,1]-h[0,0]).tolist())


def run(level,files,history,point_fixture,final_fixture):
    path=[read(p) for p in files]
    (orig,orig_sha),(last,last_sha),(alt,alt_sha),(start,start_sha)=path
    d=ScalarSubshellChart(load('fine')[0],SPECS['fine'],2,level,1e-5)
    key=dict(grid='fine',charge_div=2,scalar_div=level,cut=1e-5,counts=[d.ns,d.nf,d.np],reference_x=X)
    steps=[(.05,2048),(.05,1024),(.05,2048),(.04,1024)]
    if any((z['x_end'],z['steps'])!=q or any(z[k]!=v for k,v in key.items())
           for (z,_),q in zip(path,steps,strict=True)):
        raise ArithmeticError('The four inputs do not share one spatial chart or the required clocks')
    receipt=[z['state_sha256_uncompressed'] for z in final_fixture['full_state_panel_receipts'] if z['level']==level]
    if (final_fixture['step_counts']!=[1024,2048,1024,2048] or receipt!=[last_sha,orig_sha]
        or point_fixture['scalar_subcells']!=level or point_fixture['original_and_alternate_state_sha256']!=[orig_sha,alt_sha]
        or point_fixture['original_coarse_recovery_state_sha256']!=start_sha
        or alt['actual_prior_state_sha256']!=start_sha or history['inputs_sha256'][2]!=orig_sha
        or history['scalar_div']!=level or history['alternate_path']!=[1024,2048]
        or history['archive_prior_coarse_trace_defect']>5e-13):
        raise ArithmeticError('Input bytes do not reproduce the three separately validated corners')
    Z=[np.asarray(z['endpoint_occupation']) for z,_ in path]
    def positive(v):
        return v.shape==(d.n,) and np.all(np.isfinite(v)) and min(v)>=0 and max(v[d.ns:])<=1
    if not all(positive(v) for v in Z):
        raise ArithmeticError('Invalid source occupation')
    def count(v):
        return fsum(d.ds.w*v[:d.ns])/S
    if abs(count(Z[-1])/rows(history['earlier_coarse_recovery'])[-1,1]-1)>4e-14:
        raise ArithmeticError('Original earlier coarse state fails its paired trace')
    trace,w=drive(d,.05,1024,x_start=.04,initial=Z[-1],stepper=step)
    checked=rows(trace)
    if not positive(w) or abs(count(w)/checked[-1,1]-1)>4e-14:
        raise ArithmeticError('Combined coarse terminal state fails its verified yield')
    states=Z[:3]+[w]
    A=[optical(d,u,.05) for u in states]
    opts=[setup(d,u,v) for u,v in zip(states,A,strict=True)]
    c=[]
    for q in (3.5e-6,4.7e-6):
        p,e=zip(*(point(d,u,v,q) for u,v in zip(states,opts,strict=True)),strict=True)
        p=np.asarray(p).reshape(2,2,3)
        if np.max(abs(p[1,0]/p[0,0]-1.-next(z['alternate_over_original_minus_one'] for z in point_fixture['rows'] if z['q']==q)))>2e-12:
            raise ArithmeticError('The earlier-only edge differs from the independent archived point')
        c.append(dict(q=q,maximum_individual_mesh_defect=max(e),**decompose(p)))
    state=dict(alt,steps=1024,endpoint_occupation=w.tolist(),final_clock_control='coarse_04_to_05')
    out=dict(level=level,intervals=[[.02,.04],[.04,.05]],order_of_each_axis=[2048,1024],
             initial_states_sha256=[orig_sha,last_sha,alt_sha,start_sha],combined_terminal_trace=trace,
             scalar_yield=decompose(np.array([count(u) for u in states]).reshape(2,2)),
             optical_fits=decompose(np.asarray(A).reshape(2,2,3)),point_coefficients=c,
             columns=['forward_number','inverse_number','affine_slope'],
             scope='Matched temporal perturbations on two adjacent intervals and two fixed q points; not a whole-band or continuum time error bound.')
    return out,state


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--level',type=int,choices=(2,4),required=True)
    p.add_argument('--states',nargs=4,type=Path,required=True,metavar=('ORIGINAL','LAST_COARSE','EARLIER_COARSE','EARLIER_X04'))
    for k in ('history','points','final-fixture','combined'):
        p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args()
    js=lambda f: json.loads(f.read_text())
    pt=next(z for z in js(a.points)['results'] if z['scalar_subcells']==a.level)
    r,z=run(a.level,a.states,js(a.history),pt,js(a.final_fixture))
    t=a.combined.with_suffix('.tmp')
    wire=json.dumps(z,indent=2,allow_nan=False)+'\n'
    r['combined_terminal_state_sha256']=hashlib.sha256(wire.encode()).hexdigest()
    t.write_text(wire)
    t.replace(a.combined)
    print(json.dumps(r,indent=2,allow_nan=False))
