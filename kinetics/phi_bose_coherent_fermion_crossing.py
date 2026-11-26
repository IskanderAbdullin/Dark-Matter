"""Prescribed coherent background: fermion mass crossing and finite-lobe Dirac evolution.
All masses/momenta in GeV. No FRW, depletion, scattering or relic calculation.
"""
import math
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.special import ellipj, ellipk
m=86e-9; M1=.1; M2=.02; y=.001; delta=M1-M2
T0=2.34865e-13; gS0=3.91
n0=.12*1.05375e-5/m*(1.973269804e-14)**3

def action(A,l):
    Ic=quad(lambda th:math.cos(th)**2*math.sqrt(1+A*A*(1+math.sin(th)**2)),0,math.pi/2,epsabs=1e-12,epsrel=1e-11)[0]
    return 6*m**3/l*4*A*A*Ic/math.pi

def diag(l,A):
    phi=math.sqrt(12)*m*A/math.sqrt(l)
    rho=6*m**4/l*(A*A+A**4)
    dphi=math.sqrt(2*rho); gap=math.hypot(delta,2*y*phi)
    ml=2*(M1*M2-y*y*phi*phi)/(M1+M2+gap)
    return dict(A=A,phi=phi,rho=rho,m_light_min=ml,
                adiabatic_mass_bound=2*y*y*phi*dphi/delta/ml**2 if ml>0 else None,
                mixing_rate_over_min_rest_gap_bound=y*dphi/delta**2)

def threshold(l):
    pc=math.sqrt(M1*M2)/y
    Ac=pc*math.sqrt(l)/(math.sqrt(12)*m)
    a=(n0/action(Ac,l))**(1/3)
    return dict(Acrit=Ac,phi_crit=pc,a_crit=a,
                T_at_g106=T0/a*(gS0/106.75)**(1/3),
                T_at_g10=T0/a*(gS0/10.75)**(1/3))

def crossing(l,fac):
    th=threshold(l); pc=th['phi_crit']; A=fac*th['Acrit']
    d=diag(l,A); ph=d['phi']
    V=lambda p:m*m*p*p/2+l*p**4/24
    speed=math.sqrt(2*(d['rho']-V(pc)))
    dmdp=-2*y*y*pc/(M1+M2); ddmdp=-2*y*y*delta**2/(M1+M2)**3
    slope=abs(dmdp)*speed
    mdotdot=ddmdp*speed**2-dmdp*(m*m*pc+l*pc**3/6)
    tpair=2*quad(lambda t:1/math.sqrt(m*m+l*ph*ph*(1+math.sin(t)**2)/12),math.asin(pc/ph),math.pi/2)[0]
    n_one_cross=4*slope**1.5/(8*math.pi**3)
    return dict(pstar=math.sqrt(slope),curvature=abs(mdotdot)/slope**1.5,
                separation=tpair*math.sqrt(slope),n_one_cross=n_one_cross,
                energy_fraction_if_return_mass_M2=n_one_cross*M2/d['rho'],
                mixing_rate_gap_at_cross=y*delta*speed/(delta**2+4*y*y*pc*pc)/(M1+M2))

def prob(l,fac,p,rtol=1e-9,yuk=y):
    ph=fac*math.sqrt(M1*M2)/y
    A=ph*math.sqrt(l)/(math.sqrt(12)*m)
    ep=A*A/(1+2*A*A); K=ellipk(ep); scale=m*math.sqrt(1+2*A*A)
    M0=np.diag([M1,M2]); P=p*np.eye(2)
    def h(phi):
        M=M0+np.array([[0,yuk*phi],[yuk*phi,0]])
        return np.block([[M,P],[P,-M]])
    ev,U=np.linalg.eigh(h(0)); neg=U[:,:2].astype(complex); pos=U[:,2:]
    def f(z,u):
        phi=ph*ellipj(z,ep)[1]
        return (-1j/scale*h(phi)@u.reshape(4,2)).ravel()
    sol=solve_ivp(f,[-K,K],neg.ravel(),rtol=rtol,atol=rtol/100,method='DOP853')
    if not sol.success: raise RuntimeError(sol.message)
    V=sol.y[:,-1].reshape(4,2); B=pos.T@V
    occ=np.sum(abs(B)**2,axis=1)
    return occ.tolist(),float(np.max(abs(V.conj().T@V-np.eye(2)))),sol.nfev

if __name__=='__main__':
    import json
    result={}
    for l,A in [(4.698884e-13,2.72085),(3.82587e-8,119.995)]:
        result[str(l)]={'T20':diag(l,A),'threshold':threshold(l),
                        'crossings':{str(f):crossing(l,f) for f in [1.001,1.01,1.1,2]}}
    l=3.82587e-8
    result['finite_lobes']={}
    for fac in [1.001,1.1,2]:
        ps=crossing(l,fac)['pstar']
        vals={str(r):prob(l,fac,r*ps) for r in [0,.25,.5,1,2]}
        assert max(v[1] for v in vals.values())<1e-8
        assert all(0<=v[0][0]<=1 and 0<=v[0][1]<=1 for v in vals.values())
        assert vals['0'][0]==[0.,0.]
        result['finite_lobes'][str(fac)]=vals
    z=prob(l,1.1,.005,yuk=0)
    assert max(z[0])<1e-25
    fine=prob(l,1.1,crossing(l,1.1)['pstar'],rtol=1e-11)
    assert np.max(np.abs(np.array(fine[0])-result['finite_lobes']['1.1']['1'][0]))<1e-8
    result['zero_coupling']=z; result['refine']=fine
    print(json.dumps(result,indent=2))
