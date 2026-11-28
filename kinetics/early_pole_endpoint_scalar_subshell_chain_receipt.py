"""Verify complete two-time-mesh scalar subshell sequences and exact endpoint states."""
import argparse
import gzip
import hashlib
import json
from math import fsum,isfinite
from pathlib import Path

import numpy as np

from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart
from early_pole_endpoint_stencil_remap_pilot import SPECS,X,load
from early_pole_resolved_parent_reference import S


EPOCHS=(.0075,.01,.02,.04,.05)


def read(p):
    raw=gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes()
    return json.loads(raw),hashlib.sha256(raw).hexdigest()


def run(reports,states):
    if len(reports)!=len(EPOCHS) or len(states)!=len(EPOCHS):
        raise ArithmeticError('A complete startup and four continuations are required')
    first=read(reports[0])[0]
    tags=('grid','charge_div','scalar_div','cut','counts','reference_x')
    opt={}
    if 'warm_floor' in first:
        tags+=('warm_floor','warm_div')
        opt={k:first[k] for k in ('warm_floor','warm_div')}
    g=first['grid']
    d=ScalarSubshellChart(load(g)[0],SPECS[g],first['charge_div'],first['scalar_div'],first['cut'],**opt)
    if first['counts']!=[d.ns,d.nf,d.np] or first['reference_x']!=X:
        raise ArithmeticError('Unexpected scalar or charge chart')
    out=[]
    for p,q,lo,hi,n in zip(reports,states,(X,*EPOCHS[:-1]),EPOCHS,(2048,512,1024,1024,1024),strict=True):
        r,rhash=read(p)
        s,shash=read(q)
        if any(r.get(t)!=first[t] or s.get(t)!=first[t] for t in tags) or (
            r.get('warm_floor')!=first.get('warm_floor') or s.get('warm_floor')!=first.get('warm_floor')):
            raise ArithmeticError('The chart does not match the startup')
        if (r['starting_x'],r['x_end'],r['coarse_steps'],s['steps'],s['x_end'])!=(lo,hi,n,2*n,hi):
            raise ArithmeticError('Missing or reordered time interval')
        y=np.asarray(s['endpoint_occupation'])
        if y.shape!=(d.n,) or not np.all(np.isfinite(y)) or min(y)<0 or max(y[d.ns:])>1:
            raise ArithmeticError('Invalid saved state')
        T=r['checkpoint_trajectories']
        if len(T)!=2 or any(len(t)!=8 for t in T):
            raise ArithmeticError('Both time meshes require all eight budget samples')
        t=np.array([[r['x'] for r in t] for t in T])
        if np.max(abs(t[0]/t[1]-1))>2e-14 or abs(t[1,-1]/hi-1)>1e-14 or min(np.diff(t).ravel())<=0:
            raise ArithmeticError('Inconsistent checkpoint clocks')
        row=[v for arm in T for v in arm]
        if any(not all(isfinite(w) for w in v.values()) or v['pair_number_error']>3e-13 or
               abs(v['continuity_residual'])>3e-10 or v['max_charged']>1 for v in row):
            raise ArithmeticError('Failed declared paired-number, continuity, or Pauli budget')
        raw=T[1][-1]['scalar_yield']
        expected=fsum(d.w[:d.ns]*y[:d.ns])/S
        drift=abs(T[0][-1]['scalar_yield']/raw-1)
        if abs(expected/raw-1)>4e-14 or abs(drift-r['relative_time_drift'])>1e-14 or drift>1e-3:
            raise ArithmeticError('Saved state or time-gate mismatch')
        if abs(r['scalar_pole_mass']/d.state(hi,y)[0]-1)>4e-14 or abs(fsum(r['native_scalar_yield'])/raw-1)>4e-14:
            raise ArithmeticError('Failed scalar gap or native count aggregation')
        out.append(dict(x=hi,scalar_yield=raw,gate_temporal_relative_drift=drift,
                        max_absolute_continuity=max(abs(v['continuity_residual']) for v in row),
                        source_report_sha256=rhash,uncompressed_state_sha256=shash))
    return dict(**{t:first[t] for t in tags},validated_chain=out,
                scope='All clocks, saved states and reported budgets checked; optional diagnostic, not spatial convergence.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reports',nargs=5,type=Path,required=True)
    p.add_argument('--states',nargs=5,type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.reports,a.states),indent=2,allow_nan=False))
