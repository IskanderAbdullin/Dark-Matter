"""Run a fixed 300-node warm chart selected by the earlier single-cell source trial."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from early_pole_endpoint_gauss_warm_chain_receipt import CLOCK, EPOCHS, read, verify_chart
from early_pole_endpoint_higher_warm_gauss_history import run_chart
from early_pole_endpoint_stencil_remap_pilot import X
from early_pole_endpoint_warm_budgeted_gauss_probe import Mixed


def build(path):
    raw=path.read_bytes()
    p=json.loads(raw)
    if ((p['training_x'],p['held_check_x'])!=(.01,.05) or len(p['assigned_and_control_parent_q'])!=2
            or len(p['rows'])!=2 or [u['x'] for u in p['rows']]!=[.01,.05]):
        raise ArithmeticError('Need exact completed single-node assigned/control source trial')
    sp=[u['options'][0] for u in p['rows']]
    j=sp[0]['upgraded_parent_index']
    if any(x['upgraded_parent_index']!=j or x['total_nodes']!=300 for x in sp):
        raise ArithmeticError('No unchanged unique trained cell selection')
    orders=np.full(22,8,dtype=int)
    orders[j]=9
    d=Mixed(orders)
    taken=d.ds.ed[:-1]>=.01
    s=np.flatnonzero(taken)[j]
    if d.ds.ed[s:s+2].tolist()!=p['assigned_and_control_parent_q'][0] or d.ns!=300 or min(d.ds.w)<=0.:
        raise ArithmeticError('Held optical window or 300-node ordering changed')
    d.initial=d.lift(d.initial)
    meta=dict(grid='fine',charge_div=2,scalar_div=4,cut=1e-5,gauss_floor=.01,reference_x=X,
              nodes=[d.ns,d.nf,d.np],scalar_cell_count=d.old_count,scalar_nodal_rule='gauss-volume-v5-single-upgrade',
              fixed_rule_source_sha256=hashlib.sha256(raw).hexdigest(),warm_orders_by_cell=orders.tolist(),scalar_node_budget=300)
    if d.initial.shape!=(d.n,):
        raise ArithmeticError('Lift failed on unique trained chart')
    return d,meta


def audit(control,reports,states):
    if not 1<=len(reports)<=len(EPOCHS) or len(states)!=len(reports):
        raise ArithmeticError('Require a consecutive actual prefix')
    n=len(reports)
    d,m=build(control)
    old=None
    z=[]
    for r,s,x,y,k in zip(reports,states,(X,*EPOCHS[:-1])[:n],EPOCHS[:n],CLOCK[:n],strict=True):
        (a,ha),(b,hb)=map(read,(r,s))
        u=verify_chart(a,b,d,x,y,k,m,d.repeats)
        if a['actual_previous_raw_sha256']!=old or b['actual_previous_raw_sha256']!=old:
            raise ArithmeticError('Actual older checkpoint was not preserved')
        z.append(dict(**u,report_sha256=ha,state_sha256=hb,previous_state_sha256=old))
        old=hb
    return dict(**m,prefix=z,scope='Actual fixed 300-node chart prefix only, no extrapolation or convergence claim.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--control',type=Path,required=True)
    mode=p.add_subparsers(dest='mode',required=True)
    a=mode.add_parser('advance')
    a.add_argument('--end',type=float,required=True)
    a.add_argument('--steps',type=int,required=True)
    a.add_argument('--restart',type=Path)
    a.add_argument('--state',type=Path,required=True)
    a.add_argument('--checkpoint-prefix',type=Path)
    a.add_argument('--checkpoint-provenance')
    a.add_argument('--checkpoint-stride',type=int,default=128)
    v=mode.add_parser('audit')
    for k in ('reports','states'):
        v.add_argument('--'+k,type=Path,nargs='+',required=True)
    t=p.parse_args()
    if t.mode=='advance':
        d,m=build(t.control)
        z,s=run_chart(d,m,t.end,t.steps,t.restart,scope='Fixed 300-node warm chart selected on separate held source, actual numerical evolution only.',
                      checkpoint_prefix=t.checkpoint_prefix,checkpoint_provenance=t.checkpoint_provenance,
                      checkpoint_stride=t.checkpoint_stride)
        q=t.state.with_suffix('.stage')
        q.write_text(json.dumps(s,indent=2,allow_nan=False)+'\n')
        q.replace(t.state)
    else:
        z=audit(t.control,t.reports,t.states)
    print(json.dumps(z,indent=2,allow_nan=False))
