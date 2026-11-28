"""Independent smooth weak-moment integration versus the projected operator."""
import json

import numpy as np
from scipy.special import roots_legendre

from early_pole_scalar_elastic_stencil_reference import Elastic

M=.4
C=7.
T=.48


def f(p):
    E=np.hypot(p,M)
    return (1+.65*np.exp(-((p-.68)/.24)**2))/np.expm1((E-M+.12)/T)


def continuum(order,no):
    z,w=roots_legendre(order)
    p=C*(1+z)/2
    vol=C*w*p*p/(4*np.pi**2)
    q,k=p[:,None],p[None,:]
    aa,bb=f(p)[:,None],f(p)[None,:]
    E,Q=np.hypot(q,M),np.hypot(k,M)
    xa=q*q/(E+M)
    xb=k*k/(Q+M)
    half=(xa+xb)*.5
    d0=xa*Q+M*xb
    norm=(vol[:,None]*vol[None,:])/(128*np.pi)
    value=0.
    u,v=roots_legendre(no)
    km=C*C/(np.hypot(C,M)+M)
    for a,b in zip(u,v,strict=True):
        delta=np.maximum(d0-q*k*a,0.)
        s=4*M*M+2*delta
        flux=np.sqrt(delta*(2*M*M+delta))/(E*Q)
        cm=np.sqrt(np.maximum(q*q+k*k+2*q*k*a,0.))*np.sqrt(delta/(2*s))
        omega=.5*b*flux*norm/s
        for c,weight in zip(u,v,strict=True):
            first=np.maximum(half+cm*c,0.)
            second=np.maximum(half-cm*c,0.)
            p3=np.sqrt(first*(first+2*M))
            p4=np.sqrt(second*(second+2*M))
            cc,dd=f(p3),f(p4)
            F=aa*bb*(1+cc+dd)-cc*dd*(1+aa+bb)
            # p^4 is intentionally not a collisional invariant.
            weak=p3**4+p4**4-q**4-k**4
            value+=float(np.sum(.5*weight*omega*F*weak*((first<=km)&(second<=km))))
    return value


def projected(n,no):
    edges=np.linspace(0,C,n+1)
    p=.75*np.diff(edges**4)/np.diff(edges**3)
    d=Elastic(p,M,edges,no)
    rate=d.flux(f(p))[0]
    return float(rate@(p**4))


if __name__=='__main__':
    exact=[continuum(n,no) for n,no in [(54,14),(78,20)]]
    check=[projected(n,no) for n,no in [(20,8),(36,12),(52,16)]]
    diff=[float(c/exact[-1]-1) for c in check]
    assert abs(exact[0]/exact[1]-1)<.004
    assert np.all(np.diff(abs(np.array(diff)))<0.) and abs(diff[-1])<.04
    print(json.dumps(dict(continuum=exact,projected=check,relative_projection_errors=diff,
                         continuum_orders=[[54,14],[78,20]],projection_orders=[[20,8],[36,12],[52,16]]),
                     indent=2,allow_nan=False))
