"""Endpoint discrete gradient for the fixed-scale scalar Hartree energy."""
import json

import numpy as np
from scipy.optimize import brentq
from scipy.special import roots_legendre


def gradient(q,n0,n1,m0,m1,lam):
    """Both endpoint masses must solve the same gap, including spectators."""
    E0,E1=np.hypot(q,m0),np.hypot(q,m1)
    A=E0+E1
    D=1+.5*lam*np.sum((n0+n1)/(2*E0*E1*A))
    du=(m1-m0)*(m1+m0)
    df=(n1-n0)/E1-n0*du/(E0*E1*A)
    B=.25*du*np.sum(df/A**2)
    return .5*A+(.25*lam*B/D)*(1/E0+1/E1)


def energy(q,n,u,ub,lam):
    return float(n@np.sqrt(q*q+u)-(u-ub)**2/(2*lam))


def gap(q,n,ub,lam):
    return brentq(lambda u:u-ub-.5*lam*(n@np.reciprocal(np.sqrt(q*q+u))),
                  ub,(.02)**2,xtol=1e-25,rtol=2e-15)


def entropy(n,V):
    f=n/V
    return V*(np.log1p(f)+f*np.log1p(1/f))


def fixture():
    rng=np.random.default_rng(445)
    lam=7.65175e-8
    cases=[]
    pts,ws=roots_legendre(64)
    for domain in ((1e-6,1e-5),(1e-6,3e-3)):
        edge=np.geomspace(*domain,33)
        q=np.cbrt((edge[1:]**3+edge[:-1]**3)/2)
        V=np.diff(edge**3)/(6*np.pi**2)
        for m in (0.,7.5e-4):
            peak=np.zeros(4)
            edge_order=True
            for trial in range(80):
                a=np.exp(rng.normal(0.,1.4,len(q)))
                a*=.017/a.sum()
                b=a*np.exp(rng.uniform(-1.2,1.2,len(q)))
                if trial%2:
                    b[-8:]=a[-8:]
                ua,ub=gap(q,a,m*m,lam),gap(q,b,m*m,lam)
                delta=b-a
                d=energy(q,b,ub,m*m,lam)-energy(q,a,ua,m*m,lam)
                Q=gradient(q,a,b,np.sqrt(ua),np.sqrt(ub),lam)
                E=.5*(np.sqrt(q*q+ua)+np.sqrt(q*q+ub))
                mid=np.sqrt(q*q+gap(q,(a+b)/2,m*m,lam))
                # Independent numerical entropy secant as line quadrature.
                n=a[:,None]+.5*delta[:,None]*(1+pts)
                ds=(.5*np.log1p(V[:,None]/n)@ws)@delta
                sx=float(np.sum(entropy(b,V)-entropy(a,V)))
                peak=np.maximum(peak,abs(np.array([d-Q@delta,d-E@delta,d-mid@delta,sx-ds])))
                edge_order &= bool(np.all(np.diff(Q)>0))
            assert peak[0]<2e-19 and peak[3]<2e-23 and edge_order
            cases.append(dict(domain=domain,bare_pole=m,tests=80,ordered=edge_order,
                              max_error=dict(discrete_gradient=peak[0],endpoint_average=peak[1],
                                             midpoint=peak[2],entropy_secant=peak[3])))
    return cases


if __name__=='__main__':
    print(json.dumps(fixture(),indent=2,allow_nan=False))
