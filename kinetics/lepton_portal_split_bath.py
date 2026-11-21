"""Conditional non-isothermal photon/electron versus neutrino χ-heat pilot.

Companion to lepton_portal_neutrino_heat.py. In anomaly-free
U(1)_Le-Lmu contact candidate, vectorlike χ exchanges with both
e−/e+ and νe,νbar_e,νmu,νbar_mu. Two active ν flavors have g=4
(one helicity per particle, confirmed chiral-vector trace).
At Tγ ≪ me the physically distinct Tν,Tγ can reverse heat signs.
Tν/Tγ=(4/11)^(1/3) below is an *illustrative instantaneous*
e± annihilation limit, not numerically solved ν decoupling/BBN.
"""
import argparse
import json
import math
from kinetics.independent_massive_chie_5d import cm_heat

TNU_OVER_TGAMMA_INSTANTANEOUS=(4/11)**(1/3)


def split_bath_heat(Tgamma_me=.1,*,Tchi_over_Tgamma=.9,
                    Tnu_over_Tgamma=TNU_OVER_TGAMMA_INSTANTANEOUS,
                    Mchi_me=.1,G_unit=.001,n=28,kinetic_cut_T=30.):
    if not(all(math.isfinite(x) and x>0 for x in
              (Tgamma_me,Tchi_over_Tgamma,Tnu_over_Tgamma,Mchi_me,
               G_unit,kinetic_cut_T)) and n>=12):
        raise ValueError("valid positive scales")
    Tg=Tgamma_me
    Tc=Tg*Tchi_over_Tgamma
    Tnu=Tg*Tnu_over_Tgamma
    M=Mchi_me
    def qheat(m,Tlep):
        pmax=math.sqrt((M+kinetic_cut_T*max(Tc,Tlep))**2-M*M)
        kmax=math.sqrt((m+kinetic_cut_T*Tlep)**2-m*m)
        return cm_heat(M,m,G_unit,Tlep,Tc,pmax=pmax,kmax=kmax,
                       nr=n,nang=n,nphi=n,gchi=4.,ge=4.)
    qe=qheat(1.,Tg)
    qnu=qheat(0.,Tnu)
    return {"Tgamma_over_me":Tg,"Tchi_over_Tgamma":Tchi_over_Tgamma,
            "Tnu_over_Tgamma":Tnu_over_Tgamma,"g_e":4,"g_nu_two_flavors":4,
            "Mchi_over_me":M,"G_me_minus2":G_unit,
            "Q_e_over_G2_me9":qe/G_unit**2,
            "Q_two_neutrino_flavors_over_G2_me9":qnu/G_unit**2,
            "Q_total_over_G2_me9":(qe+qnu)/G_unit**2,
            "Qnu_over_Qe":qnu/qe if qe else None,
            "nu_temperature_is_assumed_not_derived":True,
            "not_physical_gamma_el_or_cosmological_evolution":True}


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--Tgamma-over-me",type=float,default=.1)
    p.add_argument("--Tnu-over-Tgamma",type=float,
                   default=TNU_OVER_TGAMMA_INSTANTANEOUS)
    p.add_argument("--n",type=int,default=28)
    a=p.parse_args()
    print(json.dumps(split_bath_heat(a.Tgamma_over_me,
                                    Tnu_over_Tgamma=a.Tnu_over_Tgamma,
                                    n=a.n),indent=2))
