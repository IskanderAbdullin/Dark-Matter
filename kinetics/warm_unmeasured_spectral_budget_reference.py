"""Independent headroom-to-effective-spectrum transfer on the warm clock."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.special import kve, roots_legendre
from measured_partial_hadron_reheat_fold_reference import integrands
from warm_from_empty_recycling_envelope_reference import invariant_tail_integral
from warm_late_visible_recycling_envelope_reference import (
    COUP2, VECTOR_GAMMA, VECTOR_MASS,
)
from late_radiation_background_reference import MPL,TabulatedRadiation
from three_species_decay_onsager_reference import m,M

def rate(T,elo,ehi,n):
    C=COUP2/(48*np.pi**5)
    if ehi is None:
        return C*T**9*2*invariant_tail_integral(min(100.,max(elo/T,2*m/T)))
    E=np.linspace(elo,ehi,n+1); lo,hi=E[:-1],E[1:]
    s=hi*hi
    # ds = 2 E dE. Upper all positive factors in each invariant panel.
    D=(s-VECTOR_MASS**2)**2+(VECTOR_MASS*VECTOR_GAMMA)**2
    u=lo/T
    k0,k1=[kve(j,u)*np.exp(-u) for j in (0,1)]
    q=0.
    for mass in (m,M):
        a=(hi**6)*np.sqrt(np.maximum(0.,1-4*mass**2/s))*(1+2*mass**2/s)/D
        q+=np.diff(E)@(a*(2*T*T*k1+T*hi*k0))
    return C*q

def anchor_check():
    x,w=roots_legendre(100)
    a,b=.316,.42
    E=(a+b)/2+(b-a)/2*x
    W=(b-a)*E*w
    direct=integrands(.05,[(E**2,W,np.zeros_like(W))])[0][:,1].sum()
    bound=rate(.05,a,b,400)
    z=bound/direct
    assert 1. <= z < 1.01
    return float(z)


def run(eos,prior):
    info=json.loads(Path(prior).read_bytes())
    control=info["cases"][-1]
    q=control["comparison"]
    headroom=q["extra_energy_after_including_pions_to_touch_half_GeV4"]
    TR=.05; b_lim=500.
    S=TabulatedRadiation(eos,TR); l0=np.log(S.state(b_lim)[0]); T0=S.initial
    knots=S.sg.x[(S.sg.x>l0)&(S.sg.x<T0)]
    u=np.unique(np.r_[np.linspace(l0,T0,1001),knots])[::-1]
    lg,lr=S.sg(u),S.rg(u)
    t=(S.s0-3*u-lg)/3
    h=np.exp(.5*(np.minimum(lr[:-1],lr[1:])-np.log(2))+
                (2/3)*(lg[0]-np.maximum(lg[:-1],lg[1:]))-1e-12)
    Hph=np.pi*TR*TR/(np.sqrt(45)*MPL)
    dt=np.diff(np.exp(5*t))/(5*Hph*h)
    ea=np.exp((lg[0]-np.log(2))/3)*TR
    Tlate=ea/b_lim
    assert Tlate<.01045
    ans=[]
    for lower,upper in [(0.13,.28),(.28,.316),(.316,.42),(.42,.6),(.6,.975),(.975,2.),(2.,None)]:
        # Higher invariant thermal tails are bounded with the analytic
        # width floor. At b>=500 use the same U = E_lo/Tlate lower
        # threshold and integrate a weaker T^9 tail against Hph/b^2.
        late=COUP2/(48*np.pi**5)*2*invariant_tail_integral(min(100.,lower/Tlate))*ea**9/(4*Hph*b_lim**4)
        Rs=[]
        for n in (200,400):
            q=np.array([rate(T,lower,upper,n) for T in np.exp(u[:-1])])
            J=dt@q+late
            Rs.append(dict(E_GeV4_per_unit=float(J),constant_R_ceiling_if_only_extra_band=float(headroom/J)))
        ans.append(dict(E_lower=lower,E_upper=upper,controls=Rs))
    return dict(assumed_pion_control=control["policy"],
                comoving_energy_headroom_GeV4=headroom,
                anchor_energy_upper_over_gauss_moment=anchor_check(),
                input_sha256=[hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in (eos,prior)],
                scope="Piecewise hypothetical nonnegative effective spectral caps. "
                      "No additional hadronic rate is supplied or measured.",
                cases=ans)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos")
    p.add_argument("prior")
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.prior),indent=2,allow_nan=False))
