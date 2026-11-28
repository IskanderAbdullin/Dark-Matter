"""Recover the coarse D8 terminal state to measure cold controlling-parameter clock drift."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_startup_pilot import drive
from early_pole_cell_endpoint_transverse_audit import charge_field,continuum_moments,root
from early_pole_endpoint_capacity_fast_step import step
from early_pole_endpoint_peer_charge_shape_transport import unpack
from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_resolved_parent_reference import S


def optical(d,y,x,k=.0015):
    m=d.state(x,y)[0]
    table=[]
    for mass,e,v in ((x,d.p_edges,d.wp*y[d.iP:]),(x*d.ratio,d.d_edges,d.wd*y[d.ns:d.iP])):
        p,w,_,_,_,_=continuum_moments(e,mass,v,36)
        table.append((p,w,mass))
    _,omega,residue=root(lambda b:sum(charge_field(*a,b,36) for a in table),k)
    return m,omega,residue


def run(start,end,ref):
    g=end['grid']
    d,z,m=unpack(end,g,8)
    names=('grid','subcells_per_archived_cell','counts','reference_x')
    if (start['x_end'],start['steps'],ref['starting_x'],ref['x_end'],ref['time_refinement'])!=(.04,2048,.04,.05,1024):
        raise ArithmeticError('Need the final paired clock and its archived .04 restart')
    if any(start[k]!=end[k] for k in names) or any(ref[k]!=end[k] for k in names if k!='reference_x'):
        raise ArithmeticError('Stored reference changes the D8 chart')
    u=np.asarray(start['endpoint_occupation'])
    if u.shape!=(d.n,) or not np.all(np.isfinite(u)) or min(u)<0 or max(u[d.ns:])>1:
        raise ArithmeticError('Invalid saved restart values')
    trace,y=drive(d,.05,1024,x_start=.04,initial=u,stepper=step)
    columns=('x','scalar_yield','pair_number_error','continuity_residual','max_charged')
    a,b=[np.asarray([[v[k] for k in columns] for v in t]) for t in (trace,ref['checkpoint_trajectories'][0])]
    c=ref['checkpoint_trajectories'][1][-1]['scalar_yield']
    if a.shape!=(8,5) or a.shape!=b.shape or not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
        raise ArithmeticError('Expected complete clock traces')
    if np.max(abs(a-b))>5e-13 or abs(fsum(d.ds.w*z[:d.ns])/S/c-1)>2e-14:
        raise ArithmeticError('Recovered result does not verify the archived counts and trace')
    A,B=[np.asarray(optical(d,w,.05)) for w in (y,z)]
    return dict(grid=g,charge_subcells=8,start=.04,end=.05,steps=[1024,2048],k=.0015,
                optical_columns=('scalar_gap','photon_energy','photon_residue'),
                optical_params_coarse=A.tolist(),optical_params_fine=B.tolist(),
                coarse_fine_optical_relative=(A/B-1).tolist(),
                scalar_count_coarse_fine_relative=a[-1,1]/c-1,
                maximum_archived_trace_absolute_difference=float(np.max(abs(a-b))),
                scope='Common saved D8 .04 restart; clock sensitivity of late optical parameters only.'),dict(end,endpoint_occupation=y.tolist(),steps=1024)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('start','fine','report'):
        p.add_argument(k,type=Path)
    p.add_argument('--coarse',type=Path,required=True)
    a=p.parse_args()
    r,z=run(*map(read,(a.start,a.fine,a.report)))
    s=a.coarse.with_suffix('.tmp')
    s.write_text(json.dumps(z,indent=2,allow_nan=False)+'\n')
    s.replace(a.coarse)
    print(json.dumps(r,indent=2,allow_nan=False))
