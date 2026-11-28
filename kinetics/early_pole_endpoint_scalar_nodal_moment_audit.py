"""Compare nodal and piecewise-constant scalar moments on exactly the same states."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from early_pole_endpoint_scalar_spectral_moment_budget import LABELS,integrate,load_card
from early_pole_endpoint_scalar_shape_transport import cdf,pair


CUTS=(1e-5,.01)


def reconstructed_speed(d,state,m,n):
    q=d.ds.ed
    split=np.union1d(q,CUTS)
    a,b=split[:-1],split[1:]
    ix=np.searchsorted(q,(a+b)/2,side='right')-1
    th,wt=roots_legendre(n)
    nodes=a[:,None]+(b-a)[:,None]*(th+1)/2
    x=nodes**2
    counts=(b-a)[:,None]*wt*x*state[ix,None]/(4*np.pi**2)
    speeds=counts*x/(x+m*m)
    group=np.searchsorted(CUTS,(a+b)/2,side='right')
    return np.array([[fsum(v[group==i].ravel()) for v in (counts,speeds)] for i in range(3)])


def run(paths):
    ds=[load_card(f,n) for f,n in zip(paths,(1,2,4),strict=True)]
    m=ds[1][-1]
    out=[]
    cold=[]
    for n,(d,a,_) in zip((1,2,4),ds,strict=True):
        exact=integrate(d,a,m,48)
        ref=np.array([fsum(exact[:,i]) for i in range(len(LABELS))])
        q=d.ds.p
        E=np.hypot(q,m)
        N=d.ds.w*a[:d.ns]
        node=N[:,None]*np.stack((q**0,q,q*q,q*q/(E+m),q*q/E),axis=1)
        vals=np.array([fsum(node[:,i]) for i in range(len(LABELS))])
        if abs(vals[0]/ref[0]-1)>2e-14:
            raise ArithmeticError('The two representations must preserve the same total count')
        u,v=[reconstructed_speed(d,a,m,o) for o in (24,48)]
        err=np.max(abs(u-v)/np.maximum(abs(v),1e-105))
        if err>2e-12 or abs(fsum(v[:,0])/ref[0]-1)>2e-14:
            raise ArithmeticError('Reconstructed velocity moment failed')
        cold.append(v[0,1]/v[0,0])
        out.append(dict(scalar_div=n,nodal_over_exact_minus_one=(vals/ref-1).tolist(),
                        number_weighted_velocity_squared=fsum(v[:,1])/ref[0],
                        fraction_velocity_squared_at_q_at_least_point01=v[-1,1]/fsum(v[:,1]),
                        exact_cold_conditional_mean_velocity_squared=cold[-1],quadrature_defect=float(err)))
    # For g(q)=q²/(m²+q²), g' grows on this whole cold interval.
    cut=CUTS[0]
    if cut>=m/np.sqrt(3):
        raise ArithmeticError('The cold derivative bound does not cover the interval')
    L=2*cut*m*m/(m*m+cut*cut)**2
    W=pair(*[cdf(d,f,cut) for d,f,_ in ds[1:]])['wasserstein1_q']
    delta=cold[-1]-cold[-2]
    if abs(delta)>L*W:
        raise ArithmeticError('Smooth cold weight exceeds the transport bound')
    return dict(mass=float(m),moment_columns=LABELS,rows=out,cold_cut=cut,
                exact_S4_minus_S2_cold_conditional_v_squared=delta,
                cold_lipschitz_bound_for_v_squared=L,transport_bound_for_absolute_cold_v_squared_change=L*W,
                scope='Instantaneous relativistic dispersion proxy on numerical histograms, not a cosmological velocity or pressure prediction.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('states',nargs=3,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.states),indent=2,allow_nan=False))
