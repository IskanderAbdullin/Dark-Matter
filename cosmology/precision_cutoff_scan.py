"""Precision audit of a phenomenological, fixed-abundance colder 86-eV Bose PSD.
This is NOT an actual stimulated-production mechanism.
"Half-mode" below means POWER half: P/P_CDM = 1/2, not amplitude half.
Requires classy==3.4.3.0 in the GitHub Actions runner.
"""
from __future__ import annotations
import argparse
import json
import math
import tempfile
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid
from scipy.special import zeta
from thermal_bose_class import thermal_background, params, class_spectra, first_cross

ALPHAS = (1.0, 0.5, 0.3, 0.1, 0.05, 0.03, 0.02, 0.01)
HALF_POWER_AMPLITUDE = 1/math.sqrt(2)
CRITICAL_RHO_PER_OMEGA_H2 = 2.77536627e11  # Msun / Mpc**3
# k has units 1/Mpc, NOT h/Mpc.

def psd_file(path, alpha, coarse=False):
    assert 0 < alpha <= 1
    if coarse:
        q=np.unique(np.r_[np.geomspace(3e-6,.1,900),np.linspace(.1,45.,3200)])
    else:
        q=np.unique(np.r_[np.geomspace(1e-7,.1,1500),np.linspace(.1,50.,6200)])
    f=alpha**-3 / ((2*math.pi)**3*np.expm1(q))
    np.savetxt(path,np.column_stack((q,f)),fmt="%.17e")
    number=alpha**3*4*math.pi*trapezoid(q*q*f,q)
    rad=alpha**4*4*math.pi*trapezoid(q**3*f,q)
    return {"n_points":int(len(q)),"q_min":float(q[0]),"q_max":float(q[-1]),
      "relative_number_error":float(number/(zeta(3)/math.pi**2)-1),
      "relative_radiation_error":float(rad/(alpha*math.pi**2/30)-1)}

def at(ks, arr, k):
    return float(np.interp(math.log(k),np.log(ks),arr))

def k_half_stats(ks,T,omega_m_h2):
    power=first_cross(ks,T,HALF_POWER_AMPLITUDE)
    amplitude=first_cross(ks,T,.5)
    if power is None or amplitude is None:
        raise AssertionError("CLASS k grid does not bracket one of the cutoffs")
    density=CRITICAL_RHO_PER_OMEGA_H2*omega_m_h2
    characteristic_mass=(4*math.pi/3)*density*(math.pi/power)**3
    return {"k_power_half_1_per_Mpc":power,
      "k_amplitude_half_1_per_Mpc":amplitude,
      "characteristic_power_half_mass_Msun":characteristic_mass,
      "T_at_k1":at(ks,T,1.),
      "T_at_k5":at(ks,T,5.),
      "T_at_k10":at(ks,T,10.),
      "T_at_k20":at(ks,T,20.)}

def calculate(analytic_only=False):
    b=thermal_background()
    ks=np.geomspace(.002,35.,180)
    zs=[5.]
    out={
      "hypothesis":"Unphysical unless a separate mechanism produces it: f_alpha(p)=alpha^-3/(exp(p/(alpha*T))-1)",
      "particle_mass_eV":b["m_eV"],
      "omega_h2_analytic":b["omega_h2"],
      "Delta_Neff_thermal":b["delta_Neff_early"],
      "definition":"power_half: P_Bose/P_CDM=0.5; amplitude_half: sqrt(P_Bose/P_CDM)=0.5",
      "k_unit":"1/Mpc (NOT h/Mpc)",
      "characteristic_mass":"M=(4pi/3)rho_m*(pi/k_power_half)^3: geometric diagnostic, NOT a halo mass-function limit",
      "z":5.,
      "samples":[]}
    if not analytic_only:
        from classy import Class
        pref=params(b,cold_omega=b["omega_h2"]+1e-8)
        pref["P_k_max_h/Mpc"]=65.
        cold,background=class_spectra(Class,pref,ks,zs)
        cold_P=np.asarray(cold["5.0"])
        assert np.all(cold_P>0)
        out["CDM_reference"]=background
        out["k_1_per_Mpc"]=list(map(float,ks))
    with tempfile.TemporaryDirectory() as folder:
        for alpha in ALPHAS:
            p=Path(folder)/("fine_"+str(alpha)+".dat")
            qa=psd_file(p,alpha)
            assert abs(qa["relative_number_error"])<1e-4
            assert abs(qa["relative_radiation_error"])<1e-4
            row={"alpha":alpha,"number_ratio_expected":1,
              "radiation_ratio_expected":alpha,
              "Delta_Neff_expected":alpha*b["delta_Neff_early"],
              "occupation_at_q1":alpha**-3/math.expm1(1),
              "fine_PSD":qa}
            if not analytic_only:
                pars=params(b,psd=p)
                pars["T_ncdm"]=alpha*b["T_ratio"]
                pars["P_k_max_h/Mpc"]=65.
                warm,bgw=class_spectra(Class,pars,ks,zs)
                fine_P=np.asarray(warm["5.0"])
                assert np.all(fine_P>0)
                Tfine=np.sqrt(fine_P/cold_P)
                assert np.all(np.isfinite(Tfine))
                stats=k_half_stats(ks,Tfine,background["omega_m_h2"])
                row["CLASS_background"]=bgw
                row["fine_cutoffs"]=stats
                row["linear_transfer_amplitude_T"]=list(map(float,Tfine))
                assert abs(bgw["omega_ncdm_h2"]/b["omega_h2"]-1)<.003
                coarse=Path(folder)/("coarse_"+str(alpha)+".dat")
                qb=psd_file(coarse,alpha,coarse=True)
                cp=params(b,psd=coarse)
                cp["T_ncdm"]=alpha*b["T_ratio"]
                cp["P_k_max_h/Mpc"]=65.
                coarse_data,bc=class_spectra(Class,cp,ks,zs)
                coarse_P=np.asarray(coarse_data["5.0"])
                assert np.all(coarse_P>0)
                Tc=np.sqrt(coarse_P/cold_P)
                coarse_stats=k_half_stats(ks,Tc,background["omega_m_h2"])
                rel=np.abs(coarse_P/fine_P-1)
                reliable=ks<=stats["k_power_half_1_per_Mpc"]*1.2
                row["PSD_convergence"]={
                  "coarse_PSD":qb,
                  "coarse_omega_h2":bc["omega_ncdm_h2"],
                  "fine_coarse_power_half_k_relative_shift":abs(coarse_stats["k_power_half_1_per_Mpc"]/stats["k_power_half_1_per_Mpc"]-1),
                  "max_relative_P_error_up_to_1p2_kpower":float(max(rel[reliable])),
                  "k_amplitude_half_coarse":coarse_stats["k_amplitude_half_1_per_Mpc"],
                  "max_relative_P_error_all_k_unvalidated":float(max(rel)),
                  "reliable_k_range_power":float(stats["k_power_half_1_per_Mpc"]*1.2),
                  "P_grid_effect_at_khalf":float(at(ks,rel,stats["k_power_half_1_per_Mpc"])),
                  "max_error_k_1_per_Mpc":float(ks[np.argmax(rel)])
                  }
                row["diagnostic_checks"]={
                  "power_half_k_converges_within_2_percent":row["PSD_convergence"]["fine_coarse_power_half_k_relative_shift"]<.02,
                  "pre_cutoff_P_converges_within_4_percent":row["PSD_convergence"]["max_relative_P_error_up_to_1p2_kpower"]<.04}
                print("alpha",alpha,"k_powerhalf",stats["k_power_half_1_per_Mpc"],
                      "M_char",stats["characteristic_power_half_mass_Msun"],
                      "grid",row["PSD_convergence"]["fine_coarse_power_half_k_relative_shift"],
                      flush=True)
            out["samples"].append(row)
        if not analytic_only:
            ks_power=[x["fine_cutoffs"]["k_power_half_1_per_Mpc"] for x in out["samples"]]
            assert all(b>a for a,b in zip(ks_power,ks_power[1:])),ks_power
            validations=[r["diagnostic_checks"] for r in out["samples"]]
            out["all_precision_checks_pass"]=all(all(x.values()) for x in validations)
            out["limitations"]=[
              "No physical production Boltzmann kernel for this toy f_alpha.",
              "CLASS collisionless linear perturbations only; no Lyman-alpha likelihood.",
              "Comparing a real boson with unperturbed CDM at matched omega today; early Neff differs.",
              "Even fine/coarse CLASS PSD convergence does not establish nonlinear galaxy formation.",
              "Modes at T<<0.5 can be numerically unstable; only cutoff and pre-cutoff are validated."]
    return out

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="outputs/precision_Bose_cutoffs.json")
    p.add_argument("--analytic-only",action="store_true")
    a=p.parse_args()
    result=calculate(a.analytic_only)
    path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2),encoding="utf-8")
    print("saved",path,"samples",len(result["samples"]))
    if not a.analytic_only:
        assert result["all_precision_checks_pass"], "CLASS PSD grid convergence failed: see saved output"
