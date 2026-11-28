"""Advance the frozen-budget 255-node chart without reselecting orders during evolution."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from early_pole_endpoint_gauss_warm_chain_receipt import verify_chart,CLOCK,EPOCHS,read
from early_pole_endpoint_higher_warm_gauss_history import run_chart
from early_pole_endpoint_stencil_remap_pilot import X
from early_pole_endpoint_warm_budgeted_gauss_probe import Mixed


def build(path):
    raw=path.read_bytes()
    control=json.loads(raw)
    choices=[r for r in control['choices'] if r['maximum_total_scalar_nodes']==255]
    if len(choices)!=1 or (control['training_x'],control['held_check_x'])!=(.01,.05):
        raise ArithmeticError('Require the separately accepted 255-node frozen selection')
    p=choices[0]
    rules=p['warm_orders_by_cell']
    d=Mixed(rules)
    d.initial=d.lift(d.initial)
    if (d.ns,d.nf,d.np)!=(255,448,488) or len(d.initial)!=d.n or min(d.ds.w)<=0.:
        raise ArithmeticError('Budgeted native chart or lift changed')
    e=d.ds.ed
    if e[e>=.01].tolist()!=control['warm_edges'] or np.any(d.ds.p<=np.repeat(e[:-1],d.repeats)) or np.any(d.ds.p>=np.repeat(e[1:],d.repeats)):
        raise ArithmeticError('Budgeted Gaussian nodes do not lie in the tested support')
    meta=dict(grid='fine',charge_div=2,scalar_div=4,cut=1e-5,gauss_floor=.01,reference_x=X,
              nodes=[d.ns,d.nf,d.np],scalar_cell_count=d.old_count,
              scalar_nodal_rule='gauss-volume-v4-fixed-budget',
              budget_reference_sha256=hashlib.sha256(raw).hexdigest(),
              warm_orders_by_cell=rules,scalar_node_budget=255)
    return d,meta


def audit(path,reports,states):
    d,meta=build(path)
    k=len(reports)
    if not 1<=k<=len(EPOCHS) or len(states)!=k:
        raise ArithmeticError('Need a consecutive actual chronological prefix')
    out=[]
    older=None
    for f,g,a,b,c in zip(reports,states,(X,*EPOCHS[:-1])[:k],EPOCHS[:k],CLOCK[:k],strict=True):
        r,s=[read(z) for z in (f,g)]
        row=verify_chart(r[0],s[0],d,a,b,c,meta,d.repeats)
        if any('actual_previous_raw_sha256' not in z[0] or z[0]['actual_previous_raw_sha256']!=older for z in (r,s)):
            raise ArithmeticError('Fixed budget chain failed exact restarted fine bytes')
        out.append(dict(**row,previous_state_sha256=older,report_sha256=r[1],state_sha256=s[1]))
        older=s[1]
    return dict(**meta,complete=k==len(EPOCHS),states=out,
                scope='Saved trace/state audit of fixed per-cell orders, not a physical prediction.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--control',type=Path,required=True)
    sub=p.add_subparsers(dest='op',required=True)
    A=sub.add_parser('advance')
    A.add_argument('--x-end',type=float,default=.0075)
    A.add_argument('--steps',type=int,default=2048)
    A.add_argument('--state',type=Path,required=True)
    A.add_argument('--restart',type=Path)
    B=sub.add_parser('audit')
    B.add_argument('--states',type=Path,nargs='+',required=True)
    B.add_argument('--reports',type=Path,nargs='+',required=True)
    a=p.parse_args()
    if a.op=='audit':
        r=audit(a.control,a.reports,a.states)
    else:
        d,m=build(a.control)
        r,z=run_chart(d,m,a.x_end,a.steps,a.restart,
                      scope='Fixed budget selected from a separate held source; numerical warm trajectory only.')
        q=a.state.with_suffix('.tmp')
        q.write_text(json.dumps(z,indent=2,allow_nan=False)+'\n')
        q.replace(a.state)
    print(json.dumps(r,indent=2,allow_nan=False))
