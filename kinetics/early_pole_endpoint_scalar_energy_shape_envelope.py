"""Witness residual spectral ambiguity when every shell count and mean energy is known."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from early_pole_endpoint_scalar_spectral_moment_budget import load_card
from early_pole_endpoint_scalar_shape_transport import CUT


def volume_witness(a,b,m,n):
    u,w=roots_legendre(n)
    q=a[:,None]+(b-a)[:,None]*(u+1)/2
    unit=w[None,:]*q*q
    r=q*q/(np.hypot(q,m)+m)
    return np.sum(unit*r,axis=1)/np.sum(unit,axis=1)


def bounds(a,b,m,r,n):
    R=lambda q:q*q/(np.hypot(q,m)+m)
    A,B=R(a),R(b)
    if np.any(r<=A) or np.any(r>=B):
        raise ArithmeticError('A positive mean must lie inside its shell')
    p=np.sqrt(r*(r+2*m))
    t=(r-A)/(B-A)
    # Two explicit distributions with the same mass and mean:
    # one atom at p vs two endpoint atoms at a,b.
    low=(1-t)*(p-a)+t*(b-p)
    qx,w=roots_legendre(n)
    def sample(lo,hi):
        return lo[:,None]+(hi-lo)[:,None]*(qx+1)/2
    L=R(sample(a,p))
    H=R(sample(p,b))
    # For threshold s below the mean, F<= (B-r)/(B-s).
    # Above the mean, F>= (s-r)/(s-A).
    # The integral of the pointwise width is a W1 upper bound.
    z=(p-a)/2*(((B-r)[:,None]/(B[:,None]-L))@w)
    z+=(b-p)/2*(((r-A)[:,None]/(H-A[:,None]))@w)
    if np.any(low>z+3e-14*(b-a)) or np.any(z>(b-a)*(1+3e-14)):
        raise ArithmeticError('The constructive or general support bounds are inconsistent')
    return low,z


def evaluate(d,y,m):
    a,b=d.ds.ed[:-1],d.ds.ed[1:]
    r=volume_witness(a,b,m,48)
    if np.max(abs(volume_witness(a,b,m,24)/r-1))>2e-12:
        raise ArithmeticError('Cell energy reconstruction did not resolve')
    (u,v),(l,h)=[bounds(a,b,m,r,n) for n in (30,60)]
    if np.max(abs(v/h-1))>3e-12 or not np.array_equal(u,l):
        raise ArithmeticError('CDF width integration or constructed witness changed')
    N=d.ds.w*y[:d.ns]
    def pack(mask):
        w=N[mask]/fsum(N[mask])
        return dict(constructive_W1_witness=fsum(w*l[mask]),pointwise_CDF_envelope_W1_upper=fsum(w*h[mask]),
                    same_count_support_only_W1_upper=fsum(w*(b-a)[mask]),
                    maximum_single_cell_count_fraction=float(max(w)))
    cold=b<=CUT
    if b[cold][-1]!=CUT:
        raise ArithmeticError('The cold partition must be an original cut')
    return dict(full=pack(np.ones(len(a),dtype=bool)),cold_normalized=pack(cold))


def run(paths):
    d=[load_card(f,n) for f,n in zip(paths,(1,2,4,8),strict=True)]
    m=d[1][-1]
    return dict(mass_held=float(m),cold_q_upper=CUT,
                charts=[dict(level=n,**evaluate(c,a,m)) for (c,a,_),n in zip(d,(1,2,4,8),strict=True)],
                scope='Admissible positive within-shell distributions with fixed cell means; atom examples are not evolved physical spectra.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('state',nargs=4,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.state),indent=2,allow_nan=False))
