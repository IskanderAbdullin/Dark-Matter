"""Add the small exact continuation needed to cover the lower local clock."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_direct_expansion_tracking_continuation import follow
from fd_portal_tabulated_entropy_split_expansion import Stage
from sm_dark_entropy_exchange_reference import StandardModel


def run(path,previous,clocks):
    u,v=[z.read_bytes() for z in (previous,clocks)]
    d,e=json.loads(u),json.loads(v)
    h=hashlib.sha256(path.read_bytes()).hexdigest()
    if any(z['table_input_sha256']!=h for z in (d,e)) or d['restart_and_clock_sha256'][1]!=hashlib.sha256(v).hexdigest():
        raise ArithmeticError('Tail does not match locked EOS and clock table')
    eos=StandardModel(path)
    p,q=[PositivePortal(o) for o in ((13,36,26),(17,44,32))]
    x=np.array([.91,.93,.96])
    out=[]
    for r,K in zip(d['rows'],(Stage,CombinedExpansion),strict=True):
        tag=r['operator']
        old=np.asarray(r['points'])
        if old[-1,0]!=x[0]:
            raise ArithmeticError('Expected exact .91 terminal chart')
        a,b=[follow(K(eos,o,n),old[-1],x,m,t,maximum_step=h)
             for o,n,m,t,h in ((p,112,'BDF',4e-11,.003),(q,168,'Radau',2e-12,.0012))]
        delta=np.max(abs(a[0]-b[0]),axis=0)
        if max(delta)>3e-7 or np.max(abs(old[-1]-b[0][0]))>2e-15:
            raise ArithmeticError('Independent tail did not preserve exact restart or resolve the trajectory')
        z=np.r_[old,b[0][1:]]
        t=next(y['crossing_T'] for y in e['rows'] if y['operator']==tag and y['mode']=='fast')
        ix=np.flatnonzero((z[:-1,1]>=t)&(z[1:,1]<=t))
        if len(ix)!=1:
            raise ArithmeticError('Fast local comparison temperature must be inside the exact numerical path')
        i=ix[0]
        out.append(dict(operator=tag,points=z.tolist(),fast_clock_local_T=t,fast_clock_bracket_in_log_a=z[i:i+2,0].tolist(),
                        tail_mesh_algorithm_defects=delta.tolist(),maximum_tail_three_ledger_defects=b[1].tolist()))
    return dict(rows=out,fields=d['fields'],combined_minus_timelike=(np.asarray(out[1]['points'])-out[0]['points']).tolist(),
                table_input_sha256=h,restart_and_clock_sha256=d['restart_and_clock_sha256'],
                previous_segment_sha256=hashlib.sha256(u).hexdigest(),
                scope='Joined finite restricted FD paths; temperature labels do not imply a physical reaction or decoupling threshold.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('eos','previous','clocks'):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.previous,a.clocks),indent=2,allow_nan=False))
