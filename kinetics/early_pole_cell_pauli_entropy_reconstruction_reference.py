"""Finite-quadrature positive two-moment Pauli closure for hard volume deposits."""
import json
import numpy as np
from scipy.optimize import brentq
from scipy.special import expit,roots_legendre


def fill(N,w,e,direction):
    order=np.argsort(direction*e)
    z=np.zeros_like(w)
    left=N
    for i in order:
        z[i]=min(max(left,0.),w[i])
        left-=z[i]
    return z/w


def reconstruct(N,E,w,e):
    m=w.sum()
    assert 0.<N<m
    f0,f1=(fill(N,w,e,j) for j in (1,-1))
    lo,hi=(w*e)@f0,(w*e)@f1
    tol=m*(e.max()-e.min())*3e-12
    assert lo-tol<=E<=hi+tol,(lo,E,hi)
    if E-lo<=tol:return f0,lo,hi
    if hi-E<=tol:return f1,lo,hi
    ec=(e.min()+e.max())/2
    span=(e.max()-e.min())/2
    t=(e-ec)/span
    target=(E-ec*N)/span
    a0=np.log(N/(m-N))
    def F(b):
        def M(a):return w@expit(a+b*t)-N
        a=brentq(M,a0-abs(b)-1,a0+abs(b)+1,xtol=2e-12)
        return expit(a+b*t)
    def R(b):return (w*t)@F(b)-target
    limit=1.
    while R(-limit)>0. or R(limit)<0.:
        limit*=2
        assert limit<1e10
    beta=brentq(R,-limit,limit,xtol=2e-11)
    occ=F(beta)
    assert min(occ)>=0. and max(occ)<=1.
    assert abs(w@occ-N)<m*2e-12
    assert abs((w*e)@occ-E)<m*np.max(e)*2e-12
    return occ,lo,hi


def run():
    u,w=roots_legendre(28)
    p=1.425+.023*u
    x=np.hypot(p,.05)
    e=x+2.4e-6/x
    W=w*.023*p*p
    rows=[]
    for N in (.03,.5,.92):
        total=W.sum()*N
        lo=(W*e)@fill(total,W,e,1)
        hi=(W*e)@fill(total,W,e,-1)
        for q in (.001,.1,.5,.9,.999):
            E=lo+(hi-lo)*q
            f,_,_=reconstruct(total,E,W,e)
            rows.append(dict(fill=N,position_in_energy_interval=q,min=min(f),max=max(f),
                             number_error=W@f-total,energy_error=(W*e)@f-E))
    return dict(scope='Only discrete quadrature moment feasibility, not a dynamic operator or continuum bounds.',rows=rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
