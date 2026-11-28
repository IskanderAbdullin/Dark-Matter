"""Measure same-segment pressure time-step sensitivity from accepted exact step checkpoints."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_startup_pilot import checkpoint_signature,drive_checkpoint_digest
from early_pole_endpoint_capacity_fast_step import step
from early_pole_endpoint_gauss_warm_comparison import chart_for,gauss
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart,SPECS,load


def one(ctrl,expected,base,start,end,n,initial,report,final,checkpoint,provenance):
    (a,ha),(r,hr),(f,hf),(c,hc)=[read(x) for x in (initial,report,final,checkpoint)]
    if hf!=expected or [a['x_end'],f['x_end'],r['coarse_steps']]!=[start,end,n]:
        raise ArithmeticError('Need the exact already verified restart and two-clock result')
    if (r['actual_previous_raw_sha256']!=ha or f['actual_previous_raw_sha256']!=ha or f['steps']!=2*r['coarse_steps']
            or c['completed']!=n or c['target_steps']!=n
            or c['x_reached']!=float(start*np.exp(np.log(end/start))) or c['checks']!=r['checkpoint_trajectories'][0]
            or c['payload_sha256']!=drive_checkpoint_digest(c)):
        raise ArithmeticError('Not a complete ledger-bound 1024-step snapshot')
    _,q,_=chart_for(f,budget_control=ctrl)
    key=provenance.read_text().strip()
    t=np.linspace(0.,np.log(end/start),n+1)
    if c['signature']!=checkpoint_signature(q,np.asarray(a['nodal_occupation']),t,start,end,n,step,key):
        raise ArithmeticError('Completed snapshot does not bind exact restart state and chart')
    packet=dict(f,steps=n,nodal_occupation=c['occupation'])
    m=base[1]
    A,B=[gauss(z,base[0],m,end,expected_steps=i,budget_control=ctrl)[0] for z,i in ((packet,n),(f,2*n))]
    return dict(sha256=[ha,hr,hf,hc],provenance_hash=key,coarse=A.tolist(),fine=B.tolist(),coarse_over_fine_minus_one=(A/B-1.).tolist())


def run(control,accepted,inputs):
    ab=accepted.read_bytes()
    ref=json.loads(ab)
    if len(ref['rows'])<2:
        raise ArithmeticError('Need a completely verified consecutive prefix')
    older,r=ref['rows'][-2:]
    start,end,n=older['x'],r['x'],r['steps']//2
    if (ref['charts']!=['G8','G9','mixed300'] or r['steps']!=2*n
            or (start,end,n) not in ((.01,.02,1024),(.02,.04,1024),(.04,.05,1024))):
        raise ArithmeticError('Need both actual verified three-chart states at a supported restart clock')
    d=ScalarSubshellChart(load('fine')[0],SPECS['fine'],2,4,1e-5)
    base=(d,r['common_native_mass'])
    tags=('G9','mixed300')
    binding=(None,control)
    out=[]
    for k,c,x,p in zip(tags,binding,r['states_sha256'][1:],inputs,strict=True):
        out.append(dict(chart=k,**one(c,x,base,start,end,n,*p)))
    A,B=out
    coarse=np.asarray(B['coarse'])/A['coarse']-1.
    fine=np.asarray(B['fine'])/A['fine']-1.
    start_label,end_label=(f'{x:g}'.lstrip('0') for x in (start,end))
    return dict(x=end,steps=[n,2*n],mass=base[1],fields=['count','pressure','warm_pressure'],rows=out,
                mixed300_over_G9_by_clock=[coarse.tolist(),fine.tolist()],
                fine_minus_coarse_chart_contrast=(fine-coarse).tolist(),accepted_sha256=hashlib.sha256(ab).hexdigest(),
                scope=(f'Same-start two-clock sensitivity from exact {start_label} to {end_label}. '
                       'Not a full-history time-error or continuum estimate.'))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('control','accepted'):
        p.add_argument('--'+k,type=Path,required=True)
    for k in ('g9','mixed'):
        p.add_argument('--'+k,type=Path,nargs=5,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.control,a.accepted,(a.g9,a.mixed)),indent=2,allow_nan=False))
