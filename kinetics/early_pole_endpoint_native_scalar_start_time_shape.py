"""Isolate the startup clock for the unrefined scalar spectrum on fine charge d=2."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_cell_endpoint_collision_pilot import CellChart
from early_pole_endpoint_stencil_remap_pilot import SPECS,X,load
from early_pole_resolved_parent_reference import S


def row(d,older,newer):
    x=older['x_end']
    order=(2048,4096) if x==.0075 else ((1024,1024) if x==.01 else (2048,2048))
    fields=dict(grid='fine',subcells_per_archived_cell=2,reference_x=X,
                counts=[d.ns,d.nf,d.np],x_end=x)
    got=[]
    for r,n in zip((older,newer),order,strict=True):
        if any(r.get(k)!=a for k,a in fields.items()) or r['steps']!=n:
            raise ArithmeticError('Changed fixed chart or propagation steps')
        y=np.asarray(r['endpoint_occupation'])
        if y.shape!=(d.n,) or not np.all(np.isfinite(y)) or min(y)<0 or max(y[d.ns:])>1:
            raise ArithmeticError('Invalid occupation state')
        P=d.w[:d.ns]*y[:d.ns]
        N=fsum(P)
        if N<=0:
            raise ArithmeticError('Empty scalar inventory or underflow')
        got.append((N,P/N,d.state(x,y)[0]))
    a,b=got
    z=b[1]-a[1]
    cdf=np.cumsum(z)
    ix=np.argmax(abs(cdf))
    return dict(x=x,node_counts=[d.ns,d.nf,d.np],segment_steps=order,yields=[a[0]/S,b[0]/S],
                relative_yield_change=b[0]/a[0]-1,native_normalized_tv=fsum(abs(z))/2,
                cdf_max_abs=float(abs(cdf[ix])),cdf_max_abs_upper_q=float(d.ds.ed[ix+1]),
                relative_pole_mass_change=float(b[2]/a[2]-1))


def run(lo,hi):
    expect=(.0075,.01,.02,.04,.05)
    if not 2<=len(lo)==len(hi)<=5 or [r['x_end'] for r in lo]!=list(expect[:len(lo)]):
        raise ArithmeticError('Supply chronological matched startup and continuation states')
    d=CellChart(load('fine')[0],SPECS['fine'],2,reference_x=X)
    return dict(clock_startup_steps=[2048,4096],reports=[row(d,a,b) for a,b in zip(lo,hi,strict=True)],
                scope='Only startup time refinements compared; later each segment keeps the same step count.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base',type=Path,nargs='+',required=True)
    p.add_argument('--refined-start',type=Path,nargs='+',required=True)
    a=p.parse_args()
    print(json.dumps(run([read(z) for z in a.base],[read(z) for z in a.refined_start]),indent=2,allow_nan=False))
