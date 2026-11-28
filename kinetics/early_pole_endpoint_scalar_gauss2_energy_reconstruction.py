"""Positive two-parameter occupation reconstruction for independent scalar nodal states."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp,roots_legendre

from early_pole_endpoint_scalar_gauss2_live_fixture import ScalarGauss2
from early_pole_endpoint_scalar_energy_moment_envelope import envelopes
from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_stencil_remap_pilot import SPECS,X,load


def scaled_energy(q,a,b,m):
    ea=np.hypot(a,m)
    return ((q-a)/(b-a))*((q+a)/(b+a))*((np.hypot(b,m)+ea)/(np.hypot(q,m)+ea))


def shell_rule(a,b,m,n):
    t,w=roots_legendre(n)
    q=a+(b-a)*(t+1)/2
    weights=w*q*q*(b-a)/(4*np.pi**2)
    return q,scaled_energy(q,a,b,m),np.log(weights)


def moment(S,L,beta):
    w=L+beta*S
    z=logsumexp(w)
    f=np.exp(w-z)
    return float(f@S),float(z),f


def fit(a,b,m,mass,target):
    if not mass>0 or not 0<target<1:
        raise ArithmeticError('A positive nodal population requires an interior energy')
    q,s,L=shell_rule(a,b,m,64)
    fun=lambda beta: moment(s,L,beta)[0]-target
    z=16.
    while (fun(-z)>0 or fun(z)<0) and z<2048:
        z*=2.
    if fun(-z)>0 or fun(z)<0:
        raise ArithmeticError('Need greater quadrature resolution for a nearly extremal mean')
    beta=brentq(fun,-z,z,xtol=5e-13)
    _,normal,_=moment(s,L,beta)
    norm=np.log(mass)-normal
    qa,sa,La=shell_rule(a,b,m,128)
    avg,nb,fb=moment(sa,La,beta)
    error=max(abs(np.expm1(nb+norm-np.log(mass))),abs(avg-target))
    if error>5e-12:
        raise ArithmeticError('Independent higher quadrature did not preserve both moments')
    E=np.hypot(qa,m)
    observed=np.array([fb@(qa*qa/E),fb@((qa/E)**2),fb@(1/E)])
    ends=np.array([a,b])
    r=ends*ends/(np.hypot(ends,m)+m)
    center=r[0]+target*(r[1]-r[0])
    _,box,_=envelopes(r[:1],r[1:],np.array([center]),m)
    if np.any(observed<box[0,0]-4e-14*observed) or np.any(observed>box[1,0]+4e-14*observed):
        raise ArithmeticError('Positive representative escaped independent moment envelope')
    return dict(count=mass,scaled_mean_energy=target,tilt=float(beta),
                log_occupation_at_lower_edge=float(norm),verify_error=float(error),
                mean_pressure_numerator=float(observed[0]),
                mean_velocity_squared=float(observed[1]),
                mean_inverse_energy=float(observed[2]))


def fixture(m):
    a,b=2e-7,.003
    _,S,L=shell_rule(a,b,m,128)
    tilts=(-45.,-4.,0.,3.,29.)
    rows=[fit(a,b,m,.713,moment(S,L,t)[0]) for t in tilts]
    if max(abs(r['tilt']-t) for r,t in zip(rows,tilts,strict=True))>4e-10:
        raise ArithmeticError('Independent known tilt was not recovered')


def run(c):
    g=c['grid']
    if c['scalar_nodal_rule']!='gauss-volume-v1' or c['reference_x']!=X:
        raise ArithmeticError('Expected an explicit Gauss volume nodal state')
    d=ScalarGauss2(load(g)[0],SPECS[g],c['charge_div'],c['scalar_div'],c['cut'],gauss_floor=c['gauss_floor'])
    x=c['x_end']
    y=np.asarray(c['nodal_occupation'])
    if c['nodes']!=[d.ns,d.nf,d.np] or y.shape!=(d.n,) or not np.all(np.isfinite(y)) or (
            min(y)<0 or max(y[d.ns:])>1):
        raise ArithmeticError('State does not match the declared Gauss chart')
    m=d.state(x,y)[0]
    fixture(m)
    offset=np.r_[0,np.cumsum(d.scalar_repeats)]
    N=d.ds.w*y[:d.ns]
    out=[]
    for i,(a,b) in enumerate(zip(d.ds.ed[:-1],d.ds.ed[1:],strict=True)):
        ix=slice(offset[i],offset[i+1])
        s=fsum(N[ix])
        if not s:
            continue
        t=fsum(N[ix]*scaled_energy(d.ds.p[ix],a,b,m))/s
        out.append(dict(cell=i,q_edges=[float(a),float(b)],**fit(a,b,m,s,t)))
    if abs(fsum(z['count'] for z in out)/fsum(N)-1)>3e-14:
        raise ArithmeticError('Reconstructed total count did not close')
    return dict(grid=g,x=x,mass=float(m),nonempty_cells=len(out),max_verify_error=max(z['verify_error'] for z in out),
                max_abs_tilt=max(abs(z['tilt']) for z in out),records=out,
                scope='One positive representative per shell at fixed count and mean energy; not an inferred physical shape.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('state',type=Path)
    a=p.parse_args()
    print(json.dumps(run(read(a.state)),indent=2,allow_nan=False))
