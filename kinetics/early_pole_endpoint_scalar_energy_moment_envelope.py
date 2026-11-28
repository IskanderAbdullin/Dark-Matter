"""Positive pressure and velocity envelopes at fixed shell count and mean energy."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from early_pole_endpoint_scalar_spectral_moment_budget import load_card


LABELS=('pressure_numerator','velocity_squared','inverse_energy')


def weights(r,m):
    # r=E-m. These expressions do not subtract nearly equal energies.
    p=r*(r+2*m)/(r+m)
    v=p/(r+m)
    return np.stack((p,v,1/(r+m)),axis=-1)


def envelopes(a,b,r,m):
    if min(m,np.min(a))<0 or not m>0 or np.any(a>=b) or np.any((r<a)|(r>b)):
        raise ArithmeticError('The declared mean energy must be within each cell')
    t=(r-a)/(b-a)
    p=weights(a,m)
    q=weights(b,m)
    C=np.stack((np.minimum(p,q),np.maximum(p,q)))
    chord=(1-t[:,None])*p+t[:,None]*q
    jensen=weights(r,m)
    # The stable exact chord surplus of the convex reciprocal.
    inverse_width=(r-a)*(b-r)/((m+a)*(m+b)*(m+r))
    tri=np.stack((np.column_stack((jensen[:,0]-m*m*inverse_width,chord[:,1],jensen[:,2])),
                  np.column_stack((jensen[:,:2],jensen[:,2]+inverse_width))))
    if np.any(tri[0]>tri[1]+4e-15*tri[1]):
        raise ArithmeticError('Energy-conditioned chord and Jensen bounds disagree')
    return C,tri,inverse_width


def evaluate(d,y,m):
    a,b=d.ds.ed[:-1],d.ds.ed[1:]
    E=lambda k: np.hypot(k,m)
    K=lambda k: k*k/(E(k)+m)
    lo,hi=K(a),K(b)
    mean=[]
    val=[]
    for n in (24,48):
        u,v=roots_legendre(n)
        q=a[:,None]+(b-a)[:,None]*(u+1)/2
        w=v[None,:]*q*q
        moment=w/w.sum(axis=1)[:,None]
        mean.append(np.sum(moment*K(q),axis=1))
        val.append(np.sum(moment[:,:,None]*weights(K(q),m),axis=1))
    md=np.max(abs(mean[0]/mean[1]-1))
    vd=np.max(abs(val[0]/val[1]-1))
    if max(md,vd)>2e-12:
        raise ArithmeticError('Per-cell energy or reference moments did not resolve')
    C,R,gap_width_cell=envelopes(lo,hi,mean[-1],m)
    ref=val[-1]
    if np.any(ref<R[0]-2e-14*ref) or np.any(ref>R[1]+2e-14*ref):
        raise ArithmeticError('Constant occupation test escaped the per-cell envelopes')
    N=d.ds.w*y[:d.ns]
    group=lambda z: [fsum(z[:,i]*N) for i in range(len(LABELS))]
    z=group(ref)
    c,aug=[np.array([group(j) for j in k])/np.array(z) for k in (C,R)]
    gap_width=fsum(N*gap_width_cell)
    widths=[m*m*gap_width/z[0],float(aug[1,1]-aug[0,1]),gap_width/z[2]]
    if np.max(abs(np.array(widths)-(aug[1]-aug[0])))>8e-16:
        raise ArithmeticError('Stable widths did not match the evaluated interval endpoints')
    return dict(reference_totals=z,count_only_normalized_limit=c.tolist(),
                count_and_energy_normalized_limit=aug.tolist(),
                count_only_relative_width=(c[1]-c[0]).tolist(),
                count_and_energy_relative_width=widths,
                quadrature_relative_defect=float(max(md,vd)))


def check_fixture(m):
    # A genuinely different positive distribution with support at three
    # energies must satisfy the interval, not merely the chosen histogram.
    a=np.array([0.,1e-12,1e-4,.02,3.])
    b=np.array([1e-11,1e-9,.05,2.,22.])
    q=a[:,None]+(b-a)[:,None]*np.array([.02,.21,.93])
    w=np.array([.07,.31,.62])
    _,box,_=envelopes(a,b,q@w,m)
    z=np.sum(weights(q,m)*w[None,:,None],axis=1)
    if np.any(z<box[0]-4e-15*z) or np.any(z>box[1]+4e-15*z):
        raise ArithmeticError('Three-atom support test escaped the envelopes')


def run(paths):
    lev=(1,2,4,8)
    a=[load_card(p,n) for p,n in zip(paths,lev,strict=True)]
    m=a[1][-1]
    check_fixture(m)
    return dict(mass=float(m),columns=LABELS,
                charts=[dict(level=n,**evaluate(d,y,m)) for n,(d,y,_) in zip(lev,a,strict=True)],
                scope='Conditional finite-cell expectation bounds at a common mass, not physical or solver error bounds.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('states',nargs=4,type=Path)
    x=p.parse_args()
    print(json.dumps(run(x.states),indent=2,allow_nan=False))
