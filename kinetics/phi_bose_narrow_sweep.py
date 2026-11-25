"""Conditional dilute homogeneous real-scalar swept self-resonance diagnostic.
Units eV; requires numpy/scipy. Run as a script. No physical relic prediction.
"""
import math
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
from scipy.special import expit, ellipj, ellipk

m=86.; T0=2.34865e-4; MP=1.22089e28
n0=.12*1.05375e4/m*(1.973269804e-5)**3

def thermal(T):
    z=510998.95/T
    e=lambda x: math.hypot(x,z)
    r=4/(2*math.pi**2)*quad(lambda x:x*x*e(x)*expit(-e(x)),0,80,epsabs=1e-10)[0]
    p=4/(6*math.pi**2)*quad(lambda x:x**4/e(x)*expit(-e(x)),0,80,epsabs=1e-10)[0]
    gem=2+45/(2*math.pi**2)*(r+p)
    tnu=(gem/5.5)**(1/3) if T<2e6 else 1.
    gs=gem+5.25*tnu**3
    gr=2+30/math.pi**2*r+5.25*tnu**4
    a=T0/T*((2+5.25*2/5.5)/gs)**(1/3)
    H=1.66*math.sqrt(gr)*T*T/MP
    return a,H

def diag(T,lam):
    a,H=thermal(T)
    A2=lam*n0/(6*m**3*a**3)
    I=3*math.pi/1024*m/H*A2**4
    eps=6*A2
    Ggas=math.sqrt(3)/(6144*math.pi)*lam*m*eps**3
    Gel=math.sqrt(3)/(384*math.pi)*lam*m*eps
    return dict(T=T,A=math.sqrt(A2),I=I,q=math.sqrt(3)*m*a,gas_loss_over_H=Ggas/H,hard_elastic_over_H=Gel/H,coherent_net_loss_over_H=Ggas/(24*H))

def crossing(lam,target):
    T=math.exp(brentq(lambda l: math.log(diag(math.exp(l),lam)['I']/target),math.log(1e2),math.log(1e7)))
    return diag(T,lam)

def floquet(A):
    p=A*A/(1+2*A*A); lo=4+p
    width=3*p*p/(2*math.sqrt(1-p+p*p)+2-p); h=lo+width/2; K=ellipk(p)
    def f(z,y):
        a,b,c,d=y; v=h-6*p*ellipj(z,p)[0]**2
        return [c,d,-v*a,-v*b]
    a,b,c,d=solve_ivp(f,[0,2*K],[1,0,0,1],rtol=2e-13,atol=2e-15,method='DOP853').y[:,-1]
    tr=(a+d)/math.sqrt(a*d-b*c)
    return math.sqrt(1+2*A*A)*math.acosh(abs(tr)/2)/(2*K)

def chirp(I,L):
    g=math.sqrt(2*I/math.pi)
    w=math.sqrt(L*L-g*g); b=-1j*g/(L+w)
    a=1/math.sqrt(1-abs(b)**2); b*=a
    def f(t,y): return [-1j*t*y[0]+g*y[1],g*y[0]+1j*t*y[1]]
    y=solve_ivp(f,[-L,L],[a,b],rtol=1e-9,atol=1e-11,method='DOP853').y[:,-1]
    beta=a*y[1]-1j*g/(L+w)*a*y[0]
    return abs(beta)**2,abs(y[0])**2-abs(y[1])**2

if __name__=='__main__':
    import json
    result={str(l):{str(t):crossing(l,t) for t in [.01,1,10]} for l in [4.698884e-13,3.82587e-8]}
    result['floquet_ratio']={str(A):floquet(A)/(3*A**4/32) for A in [.03,.05,.1,.2,.3]}
    result['chirp']={str(I):{str(L):chirp(I,L) for L in [20,40,80]} for I in [.01,.1,1]}
    for l in [4.698884e-13,3.82587e-8]:
        for t in [.01,1,10]:
            d=crossing(l,t); assert abs(d['I']/t-1)<1e-10 and d['A']<.03
    assert abs(result['floquet_ratio']['0.05']-1)<.01
    for I in [.01,.1,1]:
        n, norm=result['chirp'][str(I)]['80']
        assert abs(n/math.expm1(2*I)-1)<3e-6 and abs(norm-1)<2e-7
    # Swept vacuum production matches incoherent net loss divided by 4!.
    sweep=3*math.sqrt(3)/(2*math.pi**2)*2*3*math.pi/1024
    kinetic=math.sqrt(3)/(6144*math.pi)*6**4/24
    assert math.isclose(sweep,kinetic,rel_tol=1e-14)
    print(json.dumps(result,indent=2))
