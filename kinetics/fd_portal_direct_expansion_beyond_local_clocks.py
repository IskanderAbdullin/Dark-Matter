"""Continue exact direct FD histories toward the lower local rate/H equality."""
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


def run(path,archive,clocks):
    b,c=[v.read_bytes() for v in (archive,clocks)]
    d,e=json.loads(b),json.loads(c)
    h=hashlib.sha256(path.read_bytes()).hexdigest()
    tags=['timelike','combined']
    if [r['operator'] for r in d['rows']]!=tags or any(k['table_input_sha256']!=h for k in (d,e)):
        raise ArithmeticError('The exact archived models and EOS must agree')
    eos=StandardModel(path)
    p,q=[PositivePortal(o) for o in ((13,36,26),(17,44,32))]
    grid=np.array([.55,.61,.67,.73,.79,.85,.91])
    rows=[]
    for old,kind,tag in zip(d['rows'],(Stage,CombinedExpansion),tags,strict=True):
        a=np.array(old['points'][-1])
        if a[0]!=grid[0]:
            raise ArithmeticError('Restart clock changed')
        (u,v),(s,t)=[follow(kind(eos,z,n),a,grid,m,k) for z,n,m,k in ((p,112,'BDF',7e-9),(q,168,'Radau',8e-10))]
        err=np.max(abs(u-s),axis=0)
        if max(err)>3e-7:
            raise ArithmeticError('Two conservative solutions did not resolve direct continued evolution')
        bands={z['mode']:z['crossing_T'] for z in e['rows'] if z['operator']==tag}
        if set(bands)!= {'slow','fast'}:
            raise ArithmeticError('No separately computed rate/H comparison')
        table={}
        for j,T in bands.items():
            above=np.flatnonzero(s[:,1]>T); below=np.flatnonzero(s[:,1]<T)
            table[j]=None if not len(above) or not len(below) else [float(grid[max(above)]),float(grid[min(below)])]
        rows.append(dict(operator=tag,points=s.tolist(),resolved_local_clock_temperature_brackets_in_log_a=table,
                         mesh_algorithm_absolute_defects=err.tolist(),maximum_three_ledger_defects=t.tolist()))
    return dict(rows=rows,fields=d['fields'],combined_minus_timelike=(np.asarray(rows[1]['points'])-rows[0]['points']).tolist(),
                table_input_sha256=h,restart_and_clock_sha256=[hashlib.sha256(x).hexdigest() for x in (b,c)],
                scope='Finite direct restricted FD expansion only. Clock temperatures are local diagnostics, not physical events or abundance predictions.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('eos','archive','clocks'):
        p.add_argument(k,type=Path)
    t=p.parse_args()
    print(json.dumps(run(t.eos,t.archive,t.clocks),indent=2,allow_nan=False))
