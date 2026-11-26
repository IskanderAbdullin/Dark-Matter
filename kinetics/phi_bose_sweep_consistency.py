"""Formal consistency tests for unmodified homogeneous branch; not nonlinear evolution."""
from phi_bose_narrow_sweep import *

def gate(T,l):
    d=diag(T,l); A2=d['A']**2; I=d['I']
    logN=math.log(math.expm1(2*I)) if I<300 else 2*I
    logQ=math.log(math.sqrt(3)*l/(2*math.pi**2*A2))+logN
    return dict(**d,logQ=logQ,tau_band=d['hard_elastic_over_H']*A2**2/8,tau_RD=d['hard_elastic_over_H']*2*math.log(1.5))

def back(l):
    T=math.exp(brentq(lambda v:gate(math.exp(v),l)['logQ'],math.log(1e3),math.log(1e6)))
    return gate(T,l)

def jac(T):
    a1=thermal(T*math.exp(-1e-4))[0]; a2=thermal(T*math.exp(1e-4))[0]
    return math.log(a1/a2)/2e-4

def integrated(T,l):
    def f(t):
        z=math.exp(t); return math.exp(gate(z,l)['logQ'])*jac(z)
    return quad(f,math.log(100),math.log(T),epsabs=1e-9,epsrel=1e-5)[0]

def optical(T,l):
    a0,H0=thermal(T); G0=gate(T,l)['hard_elastic_over_H']
    def f(t):
        z=math.exp(t); a,H=thermal(z); u=math.sqrt(3)*a0/a; E=math.hypot(1,u)
        return G0*(a0/a)**3*H0/H*(u/(E*(1+E)))/(math.sqrt(3)/6)*jac(z)
    return quad(f,math.log(100),math.log(T),epsabs=1e-8,epsrel=1e-5)[0]

if __name__=='__main__':
    import json
    result={}
    for l in [4.698884e-13,3.82587e-8]:
        d=back(l); T=d['T']
        T2=math.exp(brentq(lambda x:integrated(math.exp(x),l)-1,math.log(T),math.log(1.1*T)))
        d2=gate(T2,l)
        assert abs(integrated(T2,l)-1)<1e-8
        result[str(l)]={'instantaneous':dict(d,integrated=integrated(T,l),optical=optical(T,l)),
                        'integrated_unity':dict(d2,Q=math.exp(d2['logQ']),optical=optical(T2,l))}
    for name,fun in [('Q1',back),('I1',lambda l:crossing(l,1))]:
        l=math.exp(brentq(lambda x:math.log(fun(math.exp(x))['q']/1.14255e-6),math.log(1e-13),math.log(1e-7)))
        result[name+'_q_reference_lambda']=l
    print(json.dumps(result,indent=2))
