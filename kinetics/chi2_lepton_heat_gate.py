"""Conditional χ2-SM heat mandated by equal-U1'-charge real-phi Yukawa branch.

For original ONE real neutral U1' scalar, y phi bar(chi2)chi1
requires q1=q2. Any diagonal vector portal of chi1 therefore
also couples chi2; the benchmark assumes both chi species
separately FD thermal with zero chemical potential and common G.
This is NOT the actual decay-produced daughter distribution.
"""
from __future__ import annotations
import argparse
import json
from kinetics.finite_mediator_5d_adaptive_reference import full_heat

M_E_GEV=.00051099895


def audit(*,M1_me=.1,M2_me=.02,mphi_eV=86.,Tlep_me=.1,
          Tchi_over_Tlep=.9,mV_me=.3,G_me_minus2=.001,
          nr=6,nmu=6,nphi=12):
    mphi_me=mphi_eV*1e-9/M_E_GEV
    if not (M1_me>M2_me+mphi_me and M2_me>0 and mphi_eV>0
            and Tlep_me>0 and Tchi_over_Tlep>0
            and mV_me>0 and G_me_minus2>0):
        raise ValueError("physical chi1->chi2+real phi threshold and parameters")
    Tc=Tlep_me*Tchi_over_Tlep
    data={}
    for name,mlep in (("two_active_nu",0.),("e_minus_e_plus",1.)):
        q1=full_heat(M=M1_me,m=mlep,Te=Tlep_me,Tchi=Tc,
                     G=G_me_minus2,mV=mV_me,nr=nr,nmu=nmu,
                     nphi=nphi,cut_kinetic_T=14.,polar="adaptive")
        q2=full_heat(M=M2_me,m=mlep,Te=Tlep_me,Tchi=Tc,
                     G=G_me_minus2,mV=mV_me,nr=nr,nmu=nmu,
                     nphi=nphi,cut_kinetic_T=14.,polar="adaptive")
        data[name]={"Q1":q1,"Q2":q2,"Q2_over_Q1":q2/q1}
    return {
       "M1_me":M1_me,"M2_me":M2_me,"mphi_me":mphi_me,
       "mV_me":mV_me,"Tlep_me":Tlep_me,"Tchi_me":Tc,
       "numerical_G_me_minus2":G_me_minus2,"nrad":nr,"nmu":nmu,
       "equal_U1_charges_q1_equals_q2_assumed":True,
       "equal_thermal_chi_FD_mu0_assumed_not_evolved":True,
       "moments":data,
       "not_actual_chi2_abundance_or_gamma_el":True,
       "not_86eV_relic_prediction":True
    }


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--M2-me",type=float,default=.02)
    p.add_argument("--nr",type=int,default=6)
    p.add_argument("--nmu",type=int,default=6)
    a=p.parse_args()
    print(json.dumps(audit(M2_me=a.M2_me,nr=a.nr,nmu=a.nmu),indent=2))
