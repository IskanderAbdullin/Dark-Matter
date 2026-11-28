"""Test a cell transplant counterfactual against independently evolved mixed300."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_gauss_warm_comparison import chart_for,gauss
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart,SPECS,load


def run(ctrl,lock,states):
    u,v=ctrl.read_bytes(),lock.read_bytes()
    control,accepted=json.loads(u),json.loads(v)
    a=accepted['rows'][-1]
    if accepted['charts']!=['G8','G9','mixed300'] or len(states)!=3:
        raise ArithmeticError('Expected exactly three independently evolved charts')
    base=ScalarSubshellChart(load('fine')[0],SPECS['fine'],2,4,1e-5)
    e=base.ds.ed
    bounds=control['assigned_and_control_parent_q'][0]
    j=int(np.searchsorted(e,bounds[0]))
    if e[j:j+2].tolist()!=bounds:
        raise ArithmeticError('The selected native parent changed')
    vals=[];cells=[];parts=[]
    for f,d,p in zip(states,a['states_sha256'],(None,None,ctrl),strict=True):
        c,h=read(f)
        if d!=h:
            raise ArithmeticError('State did not match accepted physical history')
        x,y,z=gauss(c,base,a['common_native_mass'],a['x'],expected_steps=a['steps'],budget_control=p)
        _,b,_=chart_for(c,budget_control=p)
        parts.append(np.asarray(c['nodal_occupation'])[b.ns:])
        vals.append(x)
        cells.append(z)
    if np.max(abs(np.asarray(vals[1:])/vals[0]-1.-a['ratios_to_G8']))>2e-15:
        raise ArithmeticError('Archived count and pressure ratios do not recalculate')
    n,e0=np.asarray(vals[0][:2])
    z=cells[0].copy(); z[j]=cells[1][j]
    q=(cells[2]-z)/[n,e0]
    d=(cells[2]-cells[0])/[n,e0]
    keep=np.arange(len(q))!=j
    idx=np.argmax(abs(q[keep,1]))
    h=np.flatnonzero(keep)[idx]
    others=np.array([fsum(abs(k)) for k in q[keep].T])
    absolute=np.array([fsum(abs(k)) for k in d.T])
    if (np.max(abs(np.array([fsum(k) for k in d.T])-a['ratios_to_G8'][1][:2]))>2e-14
            or min(absolute)<=0.):
        raise ArithmeticError('Actual cell sums did not close their total ratio')
    return dict(x=a['x'],control_and_actual_comparison_sha256=[hashlib.sha256(t).hexdigest() for t in (u,v)],
                selected_parent=bounds,fields=['number','pressure'],
                actual_mixed300_minus_G8_total=[fsum(k) for k in d.T],
                selected_cell_residual_vs_G9=q[j].tolist(),off_target_l1_relative_to_G8=others.tolist(),
                off_target_fraction_of_mixed300_total_l1=(others/absolute).tolist(),
                largest_off_target_pressure_parent=dict(bounds=base.ds.ed[h:h+2].tolist(),residual=q[h].tolist()),
                max_absolute_charged_occupation_differences=[np.max(abs(t-parts[0])) for t in parts[1:]],
                scope='Counterfactual diagnostic on actual evolved states at a common native mass. Not a separately evolved transplanted trajectory, exact absence of feedback, or physical prediction.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--control',type=Path,required=True)
    p.add_argument('--accepted',type=Path,required=True)
    p.add_argument('--states',type=Path,nargs=3,required=True)
    t=p.parse_args()
    print(json.dumps(run(t.control,t.accepted,t.states),indent=2,allow_nan=False))
