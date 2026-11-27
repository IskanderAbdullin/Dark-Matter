"""Cold-shape diagnostic of the optional charge-lepton annihilation portal."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.special import kve,roots_legendre

from archived_collisionless_streaming_reference import Background
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import m
from warm_late_visible_recycling_envelope_reference import VECTOR_MASS,VECTOR_GAMMA,LEPTONS,COUP2
from warm_prompt_reheating_expanding_slow_reference import visible
from warm_prompt_reheating_two_temperature_late_reference import mb_both


def rate(T,n,tail):
    x,w=roots_legendre(n)
    z=(x+1)*np.sqrt(tail)/2
    E=2*m+T*z*z
    ds=2*E*T*z*np.sqrt(tail)*w
    s=E*E
    b=np.sqrt(1-4*m*m/s)
    a=np.array(LEPTONS)[:,None]**2/s
    final=((1+2*a)*np.sqrt(np.maximum(0.,1-4*a))).sum(axis=0)
    BW=(s-VECTOR_MASS**2)**2+(VECTOR_MASS*VECTOR_GAMMA)**2
    S=COUP2*s*(1+2*m*m/s)*final/(12*np.pi*BW*b)
    assert np.all(S>=0.)
    return 4*T/(32*np.pi**4)*np.sum(ds*(s-4*m*m)*E*S*kve(1,E/T)*np.exp(-(E-2*m)/T))


def nonrelativistic_check():
    T=1e-8
    x=m/T
    n=2*m*m*T*kve(2,x)/(2*np.pi**2)
    q=4*m*m
    z=np.array(LEPTONS)**2/q
    thresh=(COUP2*q*(1+2*m*m/q)*np.sum((1+2*z)*np.sqrt(np.maximum(0.,1-4*z)))/
            (6*np.pi*((q-VECTOR_MASS**2)**2+(VECTOR_MASS*VECTOR_GAMMA)**2)))
    numerical=rate(T,72,72)/n**2
    defect=abs(numerical/thresh-1)
    assert defect<1e-5
    return dict(T=T,threshold_sigma_v_GeVminus2=float(thresh),relative_numerical_defect=float(defect))


def run(eos,prior):
    bg=Background(eos)
    d=json.loads(prior.read_text())
    rows=[]
    for i,y in enumerate(d["arms"][-1]["trajectory"]):
        if i not in (0,6,12,18,24):continue
        Tv=y["T_visible"]
        Td=y["light_scalar_T"][0]
        peak=y["phase_factors"][0]
        r,s=visible(bg,Tv)
        H=np.sqrt((r+s*sum(y["energy_per_visible_entropy"]))/3)/MPL
        c=[rate(Td,k,t) for k,t in ((52,64.),(78,78.))]
        assert abs(c[0]/c[1]-1)<1e-8
        g=mb_both(Td,52,True)[0]
        rows.append(dict(T_visible=Tv,T_light=Td,
                         vector_to_scalar_forward_gross_rate_at_same_MB_shape=c[-1]/g,
                         vector_events_per_H_s_at_held_MB_occupation=c[-1]*peak**2/(H*s),
                         scalar_events_per_H_s_at_held_MB_occupation=g*peak**2/(H*s),
                         delta_relative_quadrature=abs(c[0]/c[1]-1)))
    return dict(independent_low_energy_limit=nonrelativistic_check(),scope="Thermal-shape diagnostic only, not the separate all-invariant bound. No final blocking, "
                      "no inverse, and finite positive near-threshold integration. Does not model nonthermal high-s tails.",
                rows=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("eos",type=Path);p.add_argument("path",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.path),indent=2,allow_nan=False))
