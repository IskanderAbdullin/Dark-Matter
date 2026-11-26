"""Conditional pi0-gamma input, with a massless initial photon enhancement."""
import argparse
import hashlib
import json
import numpy as np
from pathlib import Path
from scipy.special import kve
from timelike_portal_onsager_coverage_reference import factor,alpha,V,width
from three_species_decay_onsager_reference import m,M
from late_radiation_background_reference import MPL,TabulatedRadiation
m0=.1349768
def rate(T,emax,n):
 s=np.linspace(m0*m0,emax**2,n+1); l,h=s[:-1],s[1:]; de=h-l
 El,Eh=np.sqrt(l),np.sqrt(h)
 v=El/T
 k0,k1=[kve(j,v)*np.exp(-v) for j in (0,1)]
 d=h-m0*m0
 R=np.pi*alpha/2*d**3/h**2
 R_over_d=np.pi*alpha/2*(d/h)**2
 M1=2*T*T*k1+T*Eh*k0
 M2=(6*T**3+T*h)*k1+3*T*T*Eh*k0
 term=factor*de*h**1.5/((h-V**2)**2+(V*width)**2)
 p=-1/np.expm1(-m0/T)
 out=np.zeros(2)
 for ma in (m,M):
  W=term*(h+2*ma*ma)*np.sqrt(np.maximum(0.,1-4*ma**2/h))
  out += W@np.vstack((R*M1,p*(R*M1+4*T*R_over_d*M2))).T
 return out

def cold_tail(TR,temperature_scale,b_stop,Hph,Emax):
    # Bound the subpole window by the maximum propagator and replace
    # all phase factors by 1, including two produced fermion masses.
    s=Emax**2
    rmax=np.pi*alpha/2*(s-m0*m0)**3/s**2
    rdmax=np.pi*alpha/2*(1-m0*m0/s)**2
    scale=TR*temperature_scale
    B=-1/np.expm1(-m0/(scale/b_stop))
    coef=4*factor/(V*V-s)**2/Hph
    # ds=2 E dE. ∫ u^6[2K1+uK0]=3072,
    # and ∫ u^6[(6+u^2)K1+3uK0]=27648 on [0,infinity).
    e1=coef*3072*rmax*scale**9/(4*b_stop**4)
    e2=coef*(4*27648)*rdmax*scale**11/(6*b_stop**6)
    return np.array([e1,B*(e1+e2)])


def run(eos,comparison):
    prior=json.loads(Path(comparison).read_bytes())
    budget=prior["comoving_energy_headroom_GeV4"]
    TR=.05
    stop=500.
    S=TabulatedRadiation(eos,TR)
    L,H=np.log(S.state(stop)[0]),S.initial
    K=S.sg.x[(S.sg.x>L)&(S.sg.x<H)]
    u=np.unique(np.r_[np.linspace(L,H,951),K])[::-1]
    lg,lr=S.sg(u),S.rg(u)
    t=(S.s0-3*u-lg)/3
    h=np.exp(.5*(np.minimum(lr[:-1],lr[1:])-np.log(2))+
                 (2/3)*(lg[0]-np.maximum(lg[:-1],lg[1:]))-1e-12)
    Hph=np.pi*TR*TR/(np.sqrt(45)*MPL)
    dt=np.diff(np.exp(5*t))/(5*Hph*h)
    scale=np.exp((lg[0]-np.log(2))/3+1e-12)
    assert abs(t[-1]-np.log(stop))<1e-12
    result=[]
    for top in (.28,.42,.60):
        tail=cold_tail(TR,scale,stop,Hph,top)
        cases=[]
        for n in (250,500):
            rates=np.asarray([rate(T,top,n) for T in np.exp(u[:-1])])
            assert np.all(np.diff(rates,axis=0)<=1e-180)
            E=dt@rates+tail
            assert np.all(E>0.) and E[1]>=E[0]
            cases.append(dict(s_panels=n,energy_per_unit_form_factor_squared=E.tolist(),
                              cold_tail_fraction=(tail/E).tolist(),
                              form_factor_amplitude_budget_GeVminus1=float(np.sqrt(budget/E[1]))))
        result.append(dict(invariant_E_max=top,controls=cases))
    return dict(scope="One additional neutral pion/photon vacuum channel, conditioned "
                      "on a uniform timelike form-factor magnitude cap only in the stated band. "
                      "Bose factors use massless photons and massive neutral pion at zero chemical "
                      "potential. No evidence for any assigned form-factor cap or full photon bath.",
                headroom_GeV4=budget, input_sha256=[
                    hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in (eos,comparison)],
                windows=result)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos")
    p.add_argument("comparison")
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.comparison),indent=2,allow_nan=False))
