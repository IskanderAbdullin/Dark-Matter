"""Conservative time-of-flight for two FD and two Bose cells on one matched event."""
import json
from math import fsum
from functools import reduce

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
from scipy.special import expit,logit,xlogy

S=np.array([-1.,1.,-1.,1.])
SW=np.array([1,0,3,2])


def directional(y):
    a,b,p,s=y
    return a*(1-b)*p*(1+s),(1-a)*b*(1+p)*s


def affinity(y):
    g,l=directional(y)
    if min(g,l)==0.:
        return np.inf if g>0. else (-np.inf if l>0. else 0.)
    a,b,p,s=y
    return logit(a)-logit(b)+np.log(p)-np.log1p(p)-np.log(s)+np.log1p(s)


def net(y):
    a,b=directional(y)
    return a-b


def entropy(y,w):
    a,b,p,s=y
    fb=lambda x:-xlogy(x,x)-xlogy(1-x,1-x)
    bb=lambda x:np.log1p(x)+x*np.log1p(1/x) if x>0 else 0.
    return np.dot(w,[fb(a),fb(b),bb(p),bb(s)])


def flow(y,w,t):
    y,w=np.asarray(y,dtype=float),np.asarray(w,dtype=float)
    if (not np.all(np.isfinite(np.r_[y,w,t])) or min(w)<=0 or t<0 or min(y)<0 or max(y[:2])>1):
        raise ValueError('Fermion occupations must be in [0,1], Bose occupations nonnegative')
    a=affinity(y)
    if abs(a)<1e-14 or t==0:
        return y.copy(),0.,0.
    rev=a<0
    if rev:y,w=y[SW],w[SW]
    U=min(w[0]*y[0],w[1]*(1-y[1]),w[2]*y[2])
    d=U*S/w
    def af(x):
        value=y+d*x
        if x==1.:
            return -np.inf  # one reactant is exhausted in the chosen direction
        return affinity(value)
    r=brentq(af,0.,1.,xtol=7e-17,rtol=1e-15)
    p=[np.array([z,v]) for z,v in zip(y,d,strict=True)]
    one=np.array([1.,0.])
    mul=lambda rows:reduce(np.polynomial.polynomial.polymul,rows)/U
    W=mul((p[0],one-p[1],p[2],one+p[3]))-mul((one-p[0],p[1],one+p[2],p[3]))
    H,remain=np.polynomial.polynomial.polydiv(W,[-r,1.])
    if abs(remain[0])>4e-12*max(abs(W)):
        raise ArithmeticError('Detailed balance and quartic zero did not agree')
    # u=-log(1-x/r) removes the simple time singularity at equilibrium.
    def speed(u):
        val=-np.polynomial.polynomial.polyval(-r*np.expm1(-u),H)
        if val<=0.:
            raise ArithmeticError('Rate failed its positive time-of-flight factorization')
        return val
    def E(u):
        return quad(lambda s:1/speed(s),0.,u,epsabs=6e-13,epsrel=4e-13)[0]
    hi=min(max(1.,2*t*speed(0)),65.)
    while E(hi)<t and hi<65.:
        hi=min(2*hi,65.)
    u=hi if E(hi)<t else brentq(lambda z:E(z)-t,0.,hi,xtol=8e-16,rtol=1e-15)
    x=-r*np.expm1(-u)
    z=y+d*x
    return (z[SW] if rev else z),(-U*x if rev else U*x),abs(net(z)/U-(r-x)*speed(u))


def run():
    e=np.array([.5,.8,.9,.6])  # E_P+omega = E_b+E_s
    if abs(e@S)>2e-15:raise RuntimeError('Fixture is not an on-shell energy event')
    cases=[]; thermal=[]; semigroup=[]
    for y in ([.5,.1,.9,.12],[.2,.7,.04,1.5],[.16,.36,.24,8.],[.95,.15,1.2,.02],[.38,.2,.37,.55]):
        for w in ([1.,1.,1.,1.],[.3,1.4,.7,2.],[2.,.4,1.1,.25]):
            for t in (.005,.27,1.5):
                d=lambda _,x:S*net(x)/w
                ref=solve_ivp(d,(0.,t),y,method='DOP853',atol=3e-15,rtol=5e-13).y[:,-1]
                z,I,err=flow(y,w,t)
                E=abs(fsum((I*S*e).tolist()))
                n=max(abs(np.asarray(w)*(z-y)-I*S))
                ent=entropy(z,w)-entropy(y,w)
                if ent<-2e-13 or np.any(z<=0.) or max(z[:2])>=1:
                    raise ArithmeticError('Admissibility or discrete entropy gate')
                if abs(affinity(z))>abs(affinity(y))+2e-12:
                    raise ArithmeticError('The isolated reaction did not approach equilibrium')
                z1,i1,_=flow(y,w,.4*t)
                z2,i2,_=flow(z1,w,.6*t)
                semigroup.append(max(max(abs(z-z2)),abs(I-i1-i2)))
                cases.append([max(abs(ref-z)),E,n,err,ent])
    for T in (.4,.9,1.4):
        exact=np.r_[expit((.12-e[:2])/T),1/np.expm1(e[2:]/T)]
        z,I,_=flow(exact,[1.,.8,1.6,2.],2.)
        thermal.append(max(abs(z/exact-1)))
    x=np.array(cases)
    if max(x[:,:4].ravel())>4e-12 or max(thermal)>3e-13 or max(semigroup)>3e-13:
        raise ArithmeticError((x.max(axis=0),thermal))
    stressed=[]
    for y in ((.42,.13,.8,1.06697e16),(.2,.46,.09,1.87234e15)):
        v=np.array([1.,1.3,2.2,.9])
        y=np.array(y)
        t=.7/y[-1]
        z,I,_=flow(y,v,t)
        fn=lambda _,a:np.array([net(y+S*a[0]/v)/y[-1]])
        ref=solve_ivp(fn,[0,.7],[0.],rtol=3e-13,atol=1e-15,method='DOP853').y[0,-1]
        if abs(ref-I)>5e-13:
            raise ArithmeticError(('Large-occupation extent mismatch',I,ref))
        stressed.append(dict(scalar_baseline=y[-1],physical_extent=I,observed_float64_scalar_difference=z[-1]-y[-1],
                             independent_extent_absolute_defect=abs(ref-I)))
    boundary=[]
    for st in ((.7,0.,1.,0.),(1.,0.,.4,0.),(0.,1.,0.,.9),(.0,.0,.0,2.),(.2,.3,.9,0.),(.4,.7,0.,1.2)):
        c=np.array([.62,1.43,.71,.91])
        dt=.3
        f=lambda _,z:net(z)*S/c
        a=solve_ivp(f,(0.,dt),st,method='DOP853',rtol=3e-13,atol=2e-15).y[:,-1]
        b,I,_=flow(st,c,dt)
        w=max(abs(a-b))
        if w>3e-12 or min(b)<0 or max(b[:2])>1 or entropy(b,c)+2e-13<entropy(st,c):
            raise ArithmeticError(('Boundary occupation guard',st,w))
        boundary.append(dict(initial=list(st),final=b.tolist(),explicit_extent=I,reference_error=w))
    return dict(scope='Isolated conservative quantum quartet only; no network or fitted amplitude.',
                cases=len(x),independent_dop853_abs=x[:,0].max(),energy_event_closure=x[:,1].max(),
                individual_cell_extent_closure=x[:,2].max(),quartic_factorization_abs=x[:,3].max(),
                minimum_entropy_increase=x[:,4].min(),common_temperature_stationarity=max(thermal),
                semigroup_maximum_absolute_defect=max(semigroup),large_occupation_stress=stressed,
                boundary_stress=boundary)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
