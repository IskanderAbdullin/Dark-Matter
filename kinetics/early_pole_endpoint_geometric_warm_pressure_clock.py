"""Recover a matched geometric warm-clock comparison on an unchanged .04 restart."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_startup_pilot import drive
from early_pole_endpoint_capacity_fast_step import step
from early_pole_endpoint_scalar_spectral_moment_budget import integrate
from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack
from early_pole_resolved_parent_reference import S


def run(native,restart,fine,reports,warm):
    x=.05
    d,coarse_source=unpack(restart,warm,.04)
    same,f=unpack(fine,warm,x)
    dn,yn=unpack(native,False,x)
    m=dn.state(x,yn)[0]
    if any(not np.array_equal(getattr(d.ds,k),getattr(same.ds,k)) for k in ('ed','p','w')):
        raise ArithmeticError('Source and destination grids are not identical')
    spec=('x','scalar_yield','pair_number_error','continuity_residual')
    if (reports['starting_x'],reports['x_end'],reports['coarse_steps'])!=(.04,.05,1024):
        raise ArithmeticError('Wrong archived clock segment')
    original,traced=reports['checkpoint_trajectories']
    yield_saved=fsum(d.ds.w*f[:d.ns])/S
    if abs(yield_saved/traced[-1]['scalar_yield']-1)>4e-14:
        raise ArithmeticError('Fine state does not match the archived run')
    result,y=drive(d,x,1024,x_start=.04,initial=coarse_source,stepper=step)
    z=np.array([[v[k] for k in spec] for v in result])
    h=np.array([[v[k] for k in spec] for v in original])
    diff=np.max(abs(z-h))
    if diff>5e-13 or np.max(abs(z[:,3]))>3e-10 or np.max(abs(z[:,2]))>3e-13:
        raise ArithmeticError('Recovered coarse trace does not agree with archive')
    def observable(u):
        a=integrate(d,u,m,48)
        return np.array([fsum(a[:,0]),fsum(a[:,-1]),fsum(a[-2:,-1])])
    a,b=map(observable,(y,f))
    return dict(warm_div=fine.get('warm_div',1),x_start=.04,x_end=x,steps=[1024,2048],
                max_archived_trace_difference=diff,fixed_pressure_mass=m,
                coarse_fine_relative_change=(a/b-1).tolist(),
                coarse_fine_count_pressure_and_warm_pressure=[a.tolist(),b.tolist()],
                scope='Last time segment only; the two arms use the same stored .04 state.'),dict(fine,steps=1024,endpoint_occupation=y.tolist())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('native','start','fine','report'):
        p.add_argument(k,type=Path)
    p.add_argument('--coarse',type=Path,required=True)
    p.add_argument('--warm',action='store_true')
    a=p.parse_args()
    r,v=run(*map(read,(a.native,a.start,a.fine,a.report)),a.warm)
    temp=a.coarse.with_suffix('.tmp')
    temp.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
    temp.replace(a.coarse)
    print(json.dumps(r,indent=2,allow_nan=False))
