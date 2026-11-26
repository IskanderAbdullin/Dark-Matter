"""FD-bounded two-moment entropy reference on a finite momentum cell.

An inferred negative slope is allowed as a local shape, not interpreted as
the thermodynamic temperature of the full bath. Use kinetic energy E-m.
"""
import json
import numpy as np
from scipy.optimize import brentq
from scipy.special import expit, logit, roots_legendre


def mesh(lo,hi,mass,order):
    if not (0<=lo<hi and mass>=0):
        raise ValueError("invalid momentum cell")
    z,w=roots_legendre(order)
    q=lo+(hi-lo)*(z+1)/2
    k=q*q/(np.hypot(q,mass)+mass)
    w=w*q*q
    # This exact derivative form remains stable when q=0, or when E-m
    # is much smaller than m. A numerically collapsed kinetic scale is rejected.
    l=lo*lo/(np.hypot(lo,mass)+mass) if lo else 0.
    h=hi*hi/(np.hypot(hi,mass)+mass)
    if h<=l:
        raise ValueError("collapsed kinetic-energy coordinate")
    s=(k-l)/(h-l)
    return s,w/w.sum(),l,h


def infer(lo,hi,mass,fill,mean_s,order=96):
    """Return (alpha,beta) in f=expit(alpha-beta*s), matching moment means."""
    if not (0.<fill<1. and 0.<mean_s<1.):
        raise ValueError("strictly interior normalized moments required")
    s,w,_,_=mesh(lo,hi,mass,order)
    def inner(beta):
        base=logit(fill)
        if beta==0:
            return base,fill*float(w@s)
        a,b=sorted((base,base+beta))
        def residual(alpha):
            f=expit(alpha-beta*s)
            if fill>.5:
                return (1-fill)-w@expit(beta*s-alpha)
            return w@f-fill
        alpha=brentq(residual,a,b,xtol=2e-13)
        return alpha,float(w@(s*expit(alpha-beta*s)))
    target=fill*mean_s
    a,b=-1.,1.
    for _ in range(19):
        if inner(a)[1]>=target>=inner(b)[1]:
            break
        a*=2
        b*=2
    else:
        raise ValueError("moment lies outside resolved finite-cell FD range")
    beta=brentq(lambda d:inner(d)[1]-target,a,b,xtol=2e-12)
    alpha,_=inner(beta)
    return alpha,beta


def test():
    rows=[]
    for lo,hi,m in [(0.,1e-7,.1),(3e-4,1e-3,.02),(.01,5.,.3),(0.,10.,0.)]:
        for a,b in [(-25.,3.),(.2,-5.),(0.,0.),(14.,.8),(.5,30.),(-.5,-24.)]:
            s,w,_,_=mesh(lo,hi,m,256)
            f=expit(a-b*s)
            n,k=float(w@f),float(w@(s*f))
            x,y=infer(lo,hi,m,n,k/n)
            g=expit(x-y*s)
            result=max(abs(w@g-n)/n,abs(w@(s*g)-k)/k)
            assert result<5e-11
            rows.append(dict(qlo=lo,qhi=hi,m=m,alpha=a,beta=b,
                             relative_independent_moment_residual=float(result)))
    try:
        infer(0.,1e-7,.1,.9,.99)
    except ValueError:
        pass
    else:
        raise AssertionError("impossible high-energy 90% filled cell accepted")
    return dict(cases=rows,max_relative_independent_moment_error=max(
                row["relative_independent_moment_residual"] for row in rows),
                impossible_FD_moment_pair_rejected=True)
if __name__ == "__main__":
    print(json.dumps(test(),indent=2,allow_nan=False))
