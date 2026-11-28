"""Positive two-node volume quadrature of frozen scalar cells with known constant occupation."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from early_pole_endpoint_scalar_spectral_moment_budget import load_card


LABELS=('number','q','q2','q3','kinetic_energy','pressure_numerator','velocity_squared','inverse_energy')


def volume_rule(a,b):
    # Within each q shell use q=a+(b-a)u, 0<=u<=1; normalize
    # by b to keep all moments well-conditioned for very cold cells.
    if np.any(a<0) or np.any(a>=b) or not np.all(np.isfinite(a+b)):
        raise ValueError('A nonempty finite radial cell is required')
    h=b-a
    A,H=a/b,h/b
    I=np.stack([A*A/(j+1)+2*A*H/(j+2)+H*H/(j+3) for j in range(4)],axis=-1)
    _,s,t,v=(I/I[:,:1]).T
    var=t-s*s
    if np.any(var<=0):
        raise ArithmeticError('Cell quadrature requires positive volume variance')
    J=np.stack((s,np.sqrt(var),np.sqrt(var),(v-2*s*t+s**3)/var),axis=-1).reshape(-1,2,2)
    u,V=np.linalg.eigh(J)
    W=V[:,0,:]**2
    if np.any(u<=0) or np.any(u>=1) or np.any(W<=0):
        raise ArithmeticError('An internal quadrature node or weight escaped its cell')
    volume=h*(b*b+b*a+a*a)/(6*np.pi**2)
    scale=W*volume[:,None]
    defect=max(np.max(abs(np.sum(W*u**j,axis=1)/(I[:,j]/I[:,0])-1)) for j in range(4))
    if defect>3e-14:
        raise ArithmeticError('Positive rule did not recover the first four u moments')
    return a[:,None]+h[:,None]*u,scale,float(defect)


def fields(q,m):
    E=np.hypot(q,m)
    return np.stack((q**0,q,q*q,q**3,q*q/(m+E),q*q/E,(q/E)**2,1/E),axis=-1)


def totals(a,q,w,m):
    z=a[:,None,None]*w[:,:,None]*fields(q,m)
    return np.array([fsum(z[:,:,j].ravel()) for j in range(len(LABELS))])


def result(d,y,m,level):
    a,b=d.ds.ed[:-1],d.ds.ed[1:]
    f=y[:d.ns]
    q,W,defect=volume_rule(a,b)
    d2=totals(f,q,W,m)
    base=totals(f,d.ds.p[:,None],d.ds.w[:,None],m)
    high=[]
    for n in (24,48):
        t,w=roots_legendre(n)
        k=a[:,None]+(b-a)[:,None]*(t+1)/2
        weights=(b-a)[:,None]*w[None,:]*k*k/(4*np.pi**2)
        high.append(totals(f,k,weights,m))
    ref=high[-1]
    diff=np.max(abs(high[0]/ref-1))
    if diff>3e-13 or np.max(abs(d2[:4]/ref[:4]-1))>3e-14:
        raise ArithmeticError('Reconstruction reference or polynomial identity did not close')
    return dict(level=level,bins=d.ns,positive_rule_polynomial_defect=defect,
                reference_24_to_48_relative_error=float(diff),reference_totals=ref.tolist(),
                geometric_single_node_relative_bias=(base/ref-1).tolist(),
                positive_two_node_relative_bias=(d2/ref-1).tolist())


def run(x):
    a=[load_card(f,n) for f,n in zip(x,(1,2,4,8),strict=True)]
    m=a[1][-1]
    return dict(fitted_mass=float(m),quantity_labels=LABELS,
                charts=[result(d,y,m,n) for (d,y,_),n in zip(a,(1,2,4,8),strict=True)],
                scope='Frozen scalar histogram: representation only, not dynamics or physical spectra.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('states',nargs=4,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.states),indent=2,allow_nan=False))
