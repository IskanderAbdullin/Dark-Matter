"""Full five-dimensional finite-vector Pauli heat with independent adaptive polar QUAD.

Incoming p,k,mu are tensor Gauss. For every incoming state, integrate CM
scattering polar angle by independent scipy QUAD, ordinary Gauss-cosine,
or rational Gauss-transfer; recompute boosted final-state FD Pauli
factors and physical vector trace at every event. An azimuth average
completes the 5D phase space. This tests numerical quadrature, NOT
nonthermal kinetic eigenvalues, UV validity, FRW or dark relic.
"""
from __future__ import annotations
import argparse
import json
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from kinetics.finite_mediator_forward_event_reference import integrand_factory


def gauss(a,b,n):
    x,w=leggauss(n)
    return a+(b-a)*(x+1)/2,w*(b-a)/2


def full_heat(*,M=.1,m=0.,Te=1.,Tchi=.9,G=.001,mV=.1,
              nr=8,nmu=8,nphi=16,cut_kinetic_T=14.,
              polar="adaptive",nang=32):
    if not (M>0 and m>=0 and Te>0 and Tchi>0 and G>0 and mV>0
            and nr>=4 and nmu>=4 and nphi>=8 and nang>=4
            and cut_kinetic_T>0):
        raise ValueError("physical masses, couplings, incoming temperatures and grids")
    if polar not in ("adaptive","cosine","transfer"):
        raise ValueError("polar method")
    pmax=math.sqrt((M+cut_kinetic_T*Te)**2-M*M)
    kmax=math.sqrt((m+cut_kinetic_T*Te)**2-m*m)
    p,wp=gauss(0.,pmax,nr)
    k,wk=gauss(0.,kmax,nr)
    mu,wm=gauss(-1.,1.,nmu)
    accum=0.
    if polar!="adaptive":
        z,wz=leggauss(nang)
    for pi,wi in zip(p,wp):
        E=math.hypot(pi,M)
        for ki,wki in zip(k,wk):
            omega=math.hypot(ki,m)
            for mui,wmi in zip(mu,wm):
                dot=E*omega-pi*ki*mui
                s=M*M+m*m+2*dot
                vm=math.sqrt(max(0.,dot*dot-M*M*m*m))/(E*omega)
                fun,ps2=integrand_factory(M=M,m=m,p=float(pi),k=float(ki),
                    mu=float(mui),Te=Te,Tchi=Tchi,G=G,mV=mV,nphi=nphi)
                if polar=="adaptive":
                    q=quad(fun,-1.,1.,epsabs=1e-25,
                           epsrel=2e-10,limit=300)[0]
                elif polar=="cosine":
                    q=float(np.dot(wz,[fun(c) for c in z]))
                else:
                    umax=4*ps2/(mV*mV+4*ps2)
                    u=umax*(z+1)/2
                    c=1-mV*mV*u/(2*ps2*(1-u))
                    jac=mV*mV/(2*ps2*(1-u)**2)
                    q=float(np.dot(wz,[fun(ci)*j for ci,j in zip(c,jac)])*umax/2)
                accum+=wi*pi*pi*wki*ki*ki*wmi*(vm/s)*q
    # The forward-minus-reverse weak full domain covers both orientations,
    # hence 1/2 in the 512*pi^5 prefactor (gchi=ge=4 for this benchmark).
    return accum*4*4/(512*math.pi**5)


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--mV",type=float,default=.1)
    ap.add_argument("--mlepton",type=float,default=0.)
    ap.add_argument("--nr",type=int,default=8)
    ap.add_argument("--nmu",type=int,default=8)
    ap.add_argument("--nang",type=int,default=32)
    ap.add_argument("--cut-kinetic-T",type=float,default=14.)
    a=ap.parse_args()
    result={name:full_heat(mV=a.mV,m=a.mlepton,nr=a.nr,nmu=a.nmu,
                           polar=name,nang=a.nang,
                           cut_kinetic_T=a.cut_kinetic_T)
            for name in ("adaptive","cosine","transfer")}
    result["cosine_rel_to_adaptive"]=result["cosine"]/result["adaptive"]-1
    result["transfer_rel_to_adaptive"]=result["transfer"]/result["adaptive"]-1
    result["weak_thermal_5D_only"]=True
    print(json.dumps(result,indent=2))
