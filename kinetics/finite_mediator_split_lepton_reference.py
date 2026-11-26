"""Illustrative split e±/two-active-neutrino vector-mediator thermal Pauli weak heat.

Equal squared vector charges and degeneracies 4, but unequal e/ν bath
temperatures. This is NOT actual solved neutrino decoupling, mediator
constraints, nonthermal kinetic relaxation or 86-eV relic.
"""
from __future__ import annotations
import argparse
import json
from kinetics.finite_mediator_5d_adaptive_reference import full_heat

TNU_OVER_TGAMMA=(4/11)**(1/3)


def split_heat(*,mV=.1,Tgamma=1.,Tchi=.9,
               Tnu_over_Tgamma=TNU_OVER_TGAMMA,
               M=.1,G=.001,nr=8,nmu=8,nphi=16,
               polar="adaptive",nang=32):
    if not(Tnu_over_Tgamma>0 and Tgamma>0 and Tchi>0):
        raise ValueError("positive thermal scales")
    common=dict(M=M,Tchi=Tchi,G=G,mV=mV,nr=nr,nmu=nmu,nphi=nphi,
                polar=polar,nang=nang,cut_kinetic_T=14.)
    qe=full_heat(m=1.,Te=Tgamma,**common)
    qnu=full_heat(m=0.,Te=Tnu_over_Tgamma*Tgamma,**common)
    return {"mV":mV,"Tgamma":Tgamma,"Tchi":Tchi,
            "Tnu_over_Tgamma":Tnu_over_Tgamma,
            "Qe":qe,"Qtwo_active_nu":qnu,
            "Q_total":qe+qnu,"ratio_nu_e":qnu/qe if qe else None,
            "not_real_cosmological_neutrino_decoupling":True,
            "not_physical_gamma_or_dark_relic":True}


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--mV",type=float,default=.1)
    ap.add_argument("--nr",type=int,default=8)
    ap.add_argument("--nmu",type=int,default=8)
    ap.add_argument("--method",choices=["adaptive","transfer","cosine"],
                    default="adaptive")
    a=ap.parse_args()
    print(json.dumps(split_heat(mV=a.mV,nr=a.nr,nmu=a.nmu,
                                polar=a.method),indent=2))
