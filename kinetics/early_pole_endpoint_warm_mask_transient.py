"""Advance W1 and a predetermined eight-cell warm Gauss mask one short interval."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_startup_pilot import drive
from early_pole_endpoint_capacity_fast_step import step
from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack
from early_pole_endpoint_stencil_remap_pilot import load
from early_pole_endpoint_warm_live_mask import MaskedGauss2
from early_pole_endpoint_warm_pressure_source_control import cell_weights
from early_pole_endpoint_warm_source_cells import local


def norms(d,p,delta):
    z=d.ds.w*delta[:d.ns]
    return np.array([fsum(z),fsum(p*z)])


def run(card,mask):
    x=.05
    if (mask['train_x'],mask['test_x'])!=(.01,x):
        raise ArithmeticError('Selective chart must be determined at .01, before this local test')
    a,y=unpack(card,False,x)
    row=next(r for r in mask['rows'] if r['refined_native_cells']==8)
    b=MaskedGauss2(load('fine')[0],a,row['selected_q'])
    z=b.transfer_step(y)
    m,_,_,H=a.state(x,y)
    floor=float(min(a.ds.ed[a.ds.ed>=.01]))
    P=cell_weights(a.ds.ed,m,floor,64)
    Q=np.repeat(P,b.scalar_repeats)
    i=np.repeat(b.scalar_repeats==2,b.scalar_repeats)
    q=b.ds.p
    Q[i]=q[i]**2/np.hypot(q[i],m)
    if np.max(abs(norms(b,Q,z)/norms(a,P,y)-1))>5e-11:
        raise ArithmeticError('Lift changed an initial moment')
    span=5e-5
    x_end=x*np.exp(span)
    out=[]
    for d,h,p in ((a,y,P),(b,z,Q)):
        ref=fsum(local(d,h,p,x,m,H,a.ds.ed,8)[:,2])
        arms=[]
        for n in (16,32):
            trace,end=drive(d,x_end,n,x_start=x,initial=h,stepper=step)
            if (len(trace)!=8 or max(abs(i['continuity_residual']) for i in trace)>3e-13
                    or max(i['pair_number_error'] for i in trace)>3e-13):
                raise ArithmeticError('Clock-arm conservation check failed')
            v=norms(d,p,end-h)
            if np.min(v)<=0 or not np.all(np.isfinite(v)):
                raise ArithmeticError('Expected positive short scalar moment gains')
            arms.append(dict(steps=n,number_and_pressure_gains=v.tolist(),relative_pressure_gain_to_frozen_slope=float(v[1]/(span*ref)-1),
                             max_continuity_residual=max(abs(i['continuity_residual']) for i in trace),
                             max_pair_number_error=max(i['pair_number_error'] for i in trace)))
        C,F=[np.array(a['number_and_pressure_gains']) for a in arms]
        drift=np.max(abs(C/F-1))
        if drift>3e-5:
            raise ArithmeticError(('Short clock layering did not resolve the gains',drift))
        out.append(dict(scalar_nodes=d.ns,mesh_drift=drift,frozen_pressure_slope=ref,clock_arms=arms))
    u,v=[np.asarray(r['clock_arms'][-1]['number_and_pressure_gains']) for r in out]
    return dict(x_start=x,x_end=x_end,log_span=span,fixed_measurement_mass=m,complete_warm_floor=floor,
                number_and_pressure_gain_relative_change=(v/u-1).tolist(),
                initial_frozen_pressure_slope_relative_change=out[1]['frozen_pressure_slope']/out[0]['frozen_pressure_slope']-1,
                rows=out,scope='A paired local continuation at x=.05, with fixed endpoint measurement weights but evolved state and mass. No whole-history validation.')


if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('native',type=Path)
    a.add_argument('mask',type=Path)
    o=a.parse_args()
    print(json.dumps(run(read(o.native),read(o.mask)),indent=2,allow_nan=False))
