"""FD-number, pressure and heat capacity for a chemically evolving dark bath."""
import json
import numpy as np
from scipy.optimize import brentq
from scipy.special import roots_legendre,expit
from sm_dark_entropy_exchange_reference import C
from mb_number_heat_state_reference import state as mb

def evaluate(T,eta,m=.02,*,nodes=112):
    if not (T>0 and m>0):
        raise ValueError("positive mass and temperature required")
    x,w=quad(nodes)
    # Integrate through q/T=64 for modest eta and extend before evaluating
    # the highly degenerate endpoint of a density inversion bracket.
    ymax=max(64.,eta+48.)
    u=(x+1)*np.sqrt(ymax)/2
    y=u*u
    E=T*np.hypot(m/T,y)
    a=E/T-eta
    f=expit(-a)
    ww=2*T**3/np.pi**2*w*np.sqrt(ymax)*u*y*y
    wf=f*ww
    g=wf*expit(a)
    dn=g.sum()
    n=wf.sum()
    rho=C*T**4+wf@E
    P=C*T**4/3+wf@((T*y)**2/(3*E))
    en=g@E/dn
    cv=C*4*T**3+g@((E-en)/T)**2
    S=4*C*T**3/3+ww@(f*np.logaddexp(0,a)+expit(a)*np.logaddexp(0,-a))
    return dict(n=n,rho=rho,pressure=P,entropy=S,cv_fixed_number=cv,
                chemical_energy=en,eta=eta,dn_deta=dn)

from functools import lru_cache
@lru_cache
def quad(n):
    return roots_legendre(n)

def state(T,n,m=.02,*,nodes=112):
    if n<=0.: raise ValueError("nonpositive density")
    eta_mb=mb(T,n,m)["eta"]
    p_F=(1.5*np.pi**2*n)**(1/3)
    upper=max(40.,eta_mb+15,np.hypot(p_F,m)/T+20.)
    eta=brentq(lambda u:np.log(evaluate(T,u,m,nodes=nodes)["n"]/n),
               eta_mb-1e-9,upper,xtol=1e-13)
    return evaluate(T,eta,m,nodes=nodes)

def test():
    maxerr=[]
    for T in (.005,.012,.02,.07,.25,.69,1.):
        for eta in (-13.,-3.,0.,2.,5.):
            d=evaluate(T,eta)
            r=state(T,d["n"],nodes=174)
            n=d["n"]
            ep=2e-5
            a=state(T*(1+ep),n)
            b=state(T*(1-ep),n)
            cd=(a["rho"]-b["rho"])/(2*ep*T)
            q=abs(cd/d["cv_fixed_number"]-1)
            gibbs=abs((T*(d["entropy"]+n*eta)-d["pressure"])/d["rho"]-1)
            inv=abs(r["eta"]-eta)
            md=mb(T,n)
            dil=abs(md["chemical_energy"]/d["chemical_energy"]-1)
            if eta==-13.: assert dil<3e-7
            assert max(q,gibbs,inv)<3e-8
            maxerr.append((q,gibbs,inv))
    return dict(cases=len(maxerr),max_relative_cv_check=max(r[0] for r in maxerr),
                max_relative_gibbs_check=max(r[1] for r in maxerr),
                max_absolute_eta_recovery=max(r[2] for r in maxerr),
                scope="Fixed-temperature state functions for one symmetric ideal FD species and a massless scalar, not a physical relic calculation.")
if __name__=="__main__":
    print(json.dumps(test(),indent=2,allow_nan=False))
