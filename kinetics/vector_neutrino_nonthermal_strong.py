"""Arbitrary isotropic nonthermal neutrino strong C_nu(k) from Z'↔nu nubar.

Pauli neutrino/antineutrino occupations are equal 0<f(k)<1; bosonic
b_Z'(p)>0. Exact 1↔2 shell and unpolarized massive vector gV=3,
two active ν+antiν flavors gnu=4. No ν thermal closure or BGK.
Weak vector p,cosθ and strong ν k,p integrations are independent
quadrature coordinates. No ν self-interactions, other SM, χ or FRW.
"""
from __future__ import annotations
import math
import numpy as np
from numpy.polynomial.legendre import leggauss

GV=3
GNU=4

def gl(a,b,n):
    z,w=leggauss(n)
    return a+(b-a)*(z+1)/2,w*(b-a)/2

def events(f,b,E,k,p):
    f1=np.asarray(f(k),float)
    f2=np.asarray(f(E-k),float)
    bv=np.asarray(b(p),float)
    if (np.any(f1<=0)|np.any(f1>=1)|np.any(f2<=0)|
        np.any(f2>=1)|np.any(bv<=0)|
        ~np.all(np.isfinite(f1))|~np.all(np.isfinite(f2))|
        ~np.all(np.isfinite(bv))):
        raise ValueError("strict interior Fermi and positive finite Bose modes")
    F=bv*(1-f1)*(1-f2)
    R=(1+bv)*f1*f2
    return F-R,F,R

def vector_weak(m,f,b,*,pmax=12.,np_=64,nangle=64,Gamma0=1.):
    if not(m>0 and pmax>0 and np_>=8 and nangle>=8 and Gamma0>=0):
        raise ValueError
    p,wp=gl(0.,pmax,np_)
    c,wc=gl(-1.,1.,nangle)
    E=np.hypot(p,m)
    k=(E[:,None]+p[:,None]*c[None,:])/2
    J,F,R=events(f,b,E[:,None],k,p[:,None])
    rate=Gamma0*m/E*(J@wc)/2
    w=GV/(2*math.pi**2)*wp*p*p
    dNV=-float(np.dot(w,rate))
    dEV=-float(np.dot(w,E*rate))
    entropy=float(np.sum(w[:,None]*Gamma0*m/E[:,None]*
                         wc[None,:]/2*J*np.log(F/R)))
    dSV=float(np.dot(w,rate*np.log(np.asarray(b(p))/(1+np.asarray(b(p)))))
    return {"dNV":dNV,"dEV":dEV,
            "dNnu_expected":-2*dNV,"dEnu_expected":-dEV,
            "entropy_weak":entropy,"entropy_vector":dSV}

def strong_nu_mode(k,m,f,b,*,pmax=12.,np_=64,Gamma0=1.):
    if not(m>0 and k>0 and pmax>0 and np_>=8 and Gamma0>=0):
        raise ValueError
    # Kinematics k in [(sqrt(p²+m²)-p)/2,(sqrt(p²+m²)+p)/2].
    pmin=abs(k-m*m/(4*k))
    if pmin>=pmax:return 0.
    p,wp=gl(pmin,pmax,np_)
    E=np.hypot(p,m)
    J,_,_=events(f,b,E,k,p)
    # d cosθ/2 = dk/p. 2 daughters, gV=3, gnu=4.
    return float((2*GV/GNU)*Gamma0/(k*k)*np.dot(wp,p*m/E*J))

def strong_nu_moments(m,f,b,*,pmax=12.,np_=64,nk=64,Gamma0=1.):
    if not(m>0 and pmax>0 and np_>=8 and nk>=8 and Gamma0>=0):
        raise ValueError
    Emax=math.hypot(pmax,m)
    kmin=(Emax-pmax)/2
    kmax=(Emax+pmax)/2
    # Split exact support cusp k=m/2. Unsplit Gauss is falsely biased.
    k1,w1=gl(kmin,m/2,nk)
    k2,w2=gl(m/2,kmax,nk)
    ks=np.concatenate((k1,k2))
    wk=np.concatenate((w1,w2))
    C=np.array([strong_nu_mode(float(k),m,f,b,pmax=pmax,
                    np_=np_,Gamma0=Gamma0) for k in ks])
    weights=GNU/(2*math.pi**2)*wk*ks*ks
    fk=np.asarray(f(ks),float)
    return {"dNnu":float(np.dot(weights,C)),
            "dEnu":float(np.dot(weights,ks*C)),
            "dSnu":float(np.dot(weights,C*np.log((1-fk)/fk))),
            "kmin":kmin,"kmax":kmax,
            "nonthermal_neutrino_strong_source":True,
            "no_full_neutrino_time_or_SM":True}

if __name__=="__main__":
    import json
    from scipy.special import expit
    m=.3
    f=lambda k:expit(-k+.15*np.exp(-k))
    b=lambda p:.75/np.expm1(np.hypot(p,m))
    w=vector_weak(m,f,b)
    s=strong_nu_moments(m,f,b)
    print(json.dumps({"weak":w,"strong":s,
      "Nrelative":s["dNnu"]/w["dNnu_expected"]-1,
      "Erelative":s["dEnu"]/w["dEnu_expected"]-1,
      "entropy_relative":(s["dSnu"]+w["entropy_vector"])/w["entropy_weak"]-1},indent=2))
