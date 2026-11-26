"""Thermodynamic state closure for the pending nonlinear MB number/heat test."""
import json
import numpy as np
from scipy.special import kve, roots_legendre
from sm_dark_entropy_exchange_reference import C

def state(T,n,m=.02):
    if not (T>0 and n>0 and m>0):
        raise ValueError("positive state required")
    x=m/T
    r=kve(1,x)/kve(2,x)
    neq=m*m*T*(2/np.pi**2)*kve(2,x)*np.exp(-x)
    eta=np.log(n/neq)
    mean=m*r+3*T
    dmean=3-x*x*(r*r+3*r/x-1)
    rad=C*T**4
    rho_phi, pressure_phi, cv_phi=rad,rad/3,4*rad/T
    rho=rho_phi+n*mean
    p=pressure_phi+n*T
    cv=cv_phi+n*dmean
    return dict(rho=rho,pressure=p,cv_fixed_number=cv,
                chemical_energy=mean,eta=eta,entropy=(rho+p)/T-eta*n,
                log_eq_density=float(np.log(neq)))

def test():
    checks=[]
    v,w=roots_legendre(192)
    z=np.sqrt(100.)*(v+1)/2
    y=z*z
    for T in (.006,.01,.02,.1,.3,1.):
        for n_over_T3 in (.0001,.015,.3):
            n=n_over_T3*T**3
            d=state(T,n)
            m=.02
            energy=T*np.hypot(y,m/T)
            weight=np.sqrt(100.)*z*w*y*y*np.exp(-(energy-m)/T)
            avg=weight@energy/weight.sum()
            var=weight@(energy-avg)**2/weight.sum()
            dcv=n*var/T**2+4*C*T**3
            pT=state(T*(1+2e-5),n)
            mT=state(T*(1-2e-5),n)
            a=(pT["rho"]-mT["rho"])/(4e-5*T)
            err=max(abs(avg/d["chemical_energy"]-1),abs(a/dcv-1),
                    abs(d["cv_fixed_number"]/dcv-1))
            assert err<3e-8 and d["cv_fixed_number"]>0
            checks.append(dict(T=T,n_over_T3=n_over_T3,max_relative_check=err))
    return dict(checks=checks,maximum_relative_error=max(z["max_relative_check"] for z in checks),
                status="Classical dark thermodynamics, no reaction or cosmology claim")
if __name__=="__main__":
    print(json.dumps(test(),indent=2,allow_nan=False))
