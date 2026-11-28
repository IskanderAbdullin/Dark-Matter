"""Propagate the moment-only pole enclosure into the cold kinematic threshold."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from early_pole_cell_endpoint_collision_pilot import CellChart
from early_pole_endpoint_stencil_remap_pilot import SPECS,load
from early_pole_endpoint_transverse_moment_envelope import enclosure
from finite_mass_qed_timelike_pole_batched_absorption_reference import spacelike_parent_minimum


BAND=(1e-6,1.5e-5)


def energy(M,m,u,v):
    # T^0=-u, and the admissible radial angle solves
    # v sqrt(E^2-M^2) = u E + (M^2-m^2+u^2-v^2)/2.
    s=(u-v)*(u+v)
    A=M*M+s-m*m
    if not (0<u<v and A>0):
        raise ValueError('Not on the declared monotone past spacelike branch')
    t=spacelike_parent_minimum(M,m,-u,v,s)
    def f(e):
        return v*np.sqrt((e-M)*(e+M))-u*e-A/2
    b=max(2*M,M+1.)
    while f(b)<0:
        b*=2
    check=brentq(f,M,b,xtol=3e-13,rtol=5e-15)
    if abs(t/check-1)>2e-12:
        raise ArithmeticError('Radial closed form and direct root disagree')
    return t


def run(card):
    grid=card['grid']; x=card['x_end']
    if x!=.05:
        raise ValueError('Cold band already declared only at x=.05')
    a=load(grid)[0]
    n=card['subcells_per_archived_cell']
    C=CellChart(a,SPECS[grid],n,reference_x=card['reference_x'])
    y=np.asarray(card['endpoint_occupation'])
    if card['counts']!=[C.ns,C.nf,C.np] or y.shape!=(C.n,) or min(y)<0 or max(y[C.ns:])>1:
        raise ArithmeticError('Invalid evolved state')
    m=C.state(x,y)[0]
    m1=x*C.ratio
    pole=enclosure(C,x,y)
    w0,w1=pole['tightened_frequency_enclosure']
    wm=pole['two_moment_frequency']
    k=pole['k']
    q0,q1=BAND
    E0,E1=np.hypot([q0,q1],m)
    u0,u1=w0-E1,w1-E0
    v0,v1=k-q1,k+q1
    amin=x*x-m1*m1+u0*u0-v1*v1
    margin=u0/v1-q1/E0
    if not (0<u0<=u1<v0<v1 and amin>0 and margin>0):
        raise ArithmeticError('Band is outside the analytic monotonicity proof')
    thresh=lambda w,q,v: energy(x,m1,w-np.hypot(q,m),v)
    lo,hi=thresh(w0,q1,v1),thresh(w1,q1,v0)
    profile=[thresh(wm,q1,v1),thresh(wm,q1,v0)]
    freq_only=[thresh(w,q1,np.hypot(q1,k)) for w in (w0,w1)]
    # Sample all corners and interior; this cross-check does not define the envelope.
    samples=[thresh(w,q,np.sqrt(k*k+q*q-2*k*q*c))
             for w in np.linspace(w0,w1,5) for q in np.linspace(q0,q1,9) for c in np.linspace(-1,1,9)]
    if min(samples)<lo-3e-11 or max(samples)>hi+3e-11:
        raise ArithmeticError('Threshold sweep disagrees with proved extrema')
    return dict(grid=grid,subcells=n,x=x,k=k,q_band=BAND,scalar_mass=float(m),
                frequencies=[w0,wm,w1],whole_band_radial_energy_threshold_range=[lo,hi],
                fitted_pole_angle_range=profile,q_high_cosine_zero_threshold_range=freq_only,
                q_high_cosine_zero_change_from_shape=float(freq_only[1]-freq_only[0]),
                parent_archive_max_energy=float(np.hypot(C.p_edges[-1],x)),
                spacelike_slack=float(v0-u1),monotonicity_margin=float(margin),
                scope='Kinematics only, at fixed scalar mass. Not a rate or occupation envelope.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('state',nargs='+',type=Path)
    cs=p.parse_args().state
    print(json.dumps([run(json.loads(c.read_text())) for c in cs],indent=2,allow_nan=False))
