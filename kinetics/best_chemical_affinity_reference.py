"""Optional genuine independent BEST v1.3.5 finite-Bose/Pauli affinity benchmark.

Analytic target: thermal_chemical_affinity.py integrates ALL forward and
inverse on-shell quantum factors at common T with finite chemical
potentials. Unlike an empty-daughter vacuum source, THIS tests active
Bose enhancement, Pauli blocking and reverse chi2+phi->chi1.

Default mu1=-.15, mu2=.1, mu_phi=-.3: positive affinity +.05.
Reference n_phi dot = +0.0008565027253502168 and energy rates
(chi1,phi,chi2)=(-.001810198286,.000673045188,.001137153098).
This is an instantaneous rate pass, not a self-consistent frozen PSD.

Only --execute runs third-party BEST, pinned by actual Git HEAD and
requiring optional MPI/Vegas. Ordinary tests are analytical only.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import subprocess
import sys

import numpy as np
from scipy.special import expit

from best_external_vacuum_reference import BEST_PIN
from thermal_chemical_affinity import (
    integrated_thermal_source,
    independent_physical_energy_reference,
    thermal_1to2_spectrum,
)


def affinity_protocol(mu1=-.15, mu2=.1, mu_phi=-.3,
                      grid_n=56, vegas_neval=12000):
    if (not all(np.isfinite(x) for x in (mu1,mu2,mu_phi)) or
            mu_phi>=0 or grid_n<10 or vegas_neval<1000):
        raise ValueError("Finite chemical potentials, negative Bose mu, valid grid and samples required")
    if max(abs(mu1),abs(mu2),abs(mu_phi))>5:
        raise ValueError("Common-T physical energy reference requires |mu|<=5")
    source=integrated_thermal_source(mu1=mu1,mu2=mu2,mu_phi=mu_phi)
    energy=independent_physical_energy_reference(
        mu1=mu1,mu2=mu2,mu_phi=mu_phi)
    affinity=mu1-mu2-mu_phi
    return {
        "status":"ANALYTIC PROTOCOL; independent BEST NOT EXECUTED",
        "best_pinned_sha":BEST_PIN,
        "T":1.,"mass_chi1":.1,"g1":4,"gphi":1,"gchi2":4,
        "rest_width":1.,
        "mu_chi1":mu1,"mu_chi2":mu2,"mu_phi":mu_phi,
        "chemical_affinity":affinity,
        "grid_n":grid_n,"vegas_neval":vegas_neval,
        "qmin":1e-4,"qmax":30.,
        "spin_summed_matrix_element":16*math.pi*.1*4,
        "exact_analytic_number_rate":source["n_phi_net_rate"],
        "exact_analytic_three_energy_rates":energy["physical_energy_derivatives"],
        "exact_analytic_entropy_production":affinity*source["n_phi_net_rate"],
        "independent_2D_energy_check":energy["energy_conservation_residual"],
        "external_best_executed":False,
        "self_consistent_frw_or_class":False,
    }


def execute_best(best_root,output,mu1=-.15,mu2=.1,mu_phi=-.3,
                 grid_n=56,vegas_neval=12000):
    p=affinity_protocol(mu1,mu2,mu_phi,grid_n,vegas_neval)
    folder=Path(best_root).resolve()
    if not (folder/"besthep.py").is_file():
        raise ValueError("Independent BEST besthep.py not found")
    try:
        actual=subprocess.run(
            ["git","-C",str(folder),"rev-parse","HEAD"],
            check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
            timeout=10,text=True).stdout.strip()
    except (OSError,subprocess.SubprocessError) as error:
        raise ValueError("Cannot verify independent BEST Git provenance") from error
    if actual!=BEST_PIN:
        raise ValueError("BEST Git SHA differs from pinned independent source: "+actual)
    sys.path.insert(0,str(folder))
    from besthep import BEST,grid_weights  # optional external solver
    best=BEST(q_min=p["qmin"],q_max=p["qmax"],n_grid=grid_n)
    m=p["mass_chi1"]
    best.initialize_species(
        "chi1",lambda q:expit(-(np.hypot(q,m)-mu1)),
        mass=m,dof=4,stat="fermion")
    best.initialize_species(
        "phi",lambda q:1/np.expm1(np.asarray(q)-mu_phi),
        mass=0.,dof=1,stat="boson")
    best.initialize_species(
        "chi2",lambda q:expit(-(np.asarray(q)-mu2)),
        mass=0.,dof=4,stat="fermion")
    best.add_process(
        "dec",["chi1"],["phi","chi2"],
        lambda mom:np.full(mom.shape[2],p["spin_summed_matrix_element"]),
        neval=vegas_neval,nitn=2)
    raw=best.compute_rates(["dec"])
    if best.world_rank!=0:
        return None
    rates={}
    shapes={}
    entropy_rate=0.
    for species,dof,mass,mu in (
        ("chi1",4,m,mu1),("phi",1,0.,mu_phi),("chi2",4,0.,mu2)):
        q=np.asarray(best.r_grids[species],float)
        net=np.asarray(raw[species],float)
        pref=dof/(2*np.pi**2)*grid_weights(q)*q**2
        E=np.hypot(q,mass)
        number_rate=float(np.dot(pref,net))
        energy_rate=float(np.dot(pref*E,net))
        entropy_rate+=float(np.dot(pref*(E-mu),net))
        gain=np.asarray(best.gain_rates[species],float)
        loss=np.asarray(best.loss_rates[species],float)
        gross=float(np.dot(pref,gain+loss))
        rates[species]={
            "number_rate":number_rate,"energy_rate":energy_rate,
            "gross_gain_plus_loss_number_rate":gross,
            "net_over_gross":number_rate/max(gross,1e-300),
        }
        if species=="phi":
            ref=thermal_1to2_spectrum(q,mu1=mu1,mu2=mu2,mu_phi=mu_phi)
            observed=q**2/(2*np.pi**2)*net
            peak=float(np.max(abs(ref)))
            if peak<=0:raise ValueError("Equilibrium zero source is unsuitable for affinity shape comparison")
            support=abs(ref)>.001*peak
            if support.sum()<3:raise ValueError("Insufficient q support for physical affinity profile")
            shapes={
                "source_q":q.tolist(),
                "best_signed_dR_dq":observed.tolist(),
                "analytic_signed_dR_dq":ref.tolist(),
                "robust_q_count":int(support.sum()),
                "source_shape_peak_normalized_rms":
                    float(np.sqrt(np.mean(((observed[support]-ref[support])/peak)**2))),
                "support_max_relative_difference":
                    float(np.max(abs((observed[support]-ref[support])/ref[support])))
            }
    reference=p["exact_analytic_number_rate"]
    refE=p["exact_analytic_three_energy_rates"]
    dE={name:float(rates[name]["energy_rate"]-refE[name])
        for name in ("chi1","phi","chi2")}
    result={
        "independent_best_actually_executed":True,
        "best_actual_pinned_sha":actual,
        "protocol":p,
        "all_three_species_rates":rates,
        "phi_differential_source":shapes,
        "analytic_n_phi_reference":reference,
        "relative_phi_number_source_error":rates["phi"]["number_rate"]/reference-1,
        "relative_all_three_energy_rate_errors":{
            k:dE[k]/max(abs(refE[k]),1e-300) for k in dE},
        "total_number_charge_residuals":{
            "chi1_plus_phi":rates["chi1"]["number_rate"]+rates["phi"]["number_rate"],
            "chi1_plus_chi2":rates["chi1"]["number_rate"]+rates["chi2"]["number_rate"],
        },
        "net_total_physical_energy_rate":sum(v["energy_rate"] for v in rates.values()),
        "observed_collisional_entropy_rate":entropy_rate,
        "analytic_collisional_entropy_rate":p["exact_analytic_entropy_production"],
        "external_very_first_nonequilibrium_bose_inverse_decay_only":True,
        "nonlinear_time_evolution_or_cosmological_DM_validated":False,
        "warning":"Finite common-T snapshot only; does NOT establish Bose-stimulated FRW freezeout, zero-mode, full energy closure or CLASS."
    }
    file=Path(output);file.parent.mkdir(parents=True,exist_ok=True)
    file.write_text(json.dumps(result,indent=2),encoding="utf-8")
    print(json.dumps(result,indent=2))
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--execute",action="store_true")
    parser.add_argument("--best-root",default="")
    parser.add_argument("--output",default="outputs/best_chemical_affinity.json")
    parser.add_argument("--mu1",type=float,default=-.15)
    parser.add_argument("--grid-n",type=int,default=56)
    parser.add_argument("--neval",type=int,default=12000)
    args=parser.parse_args()
    if args.execute:
        if not args.best_root:raise ValueError("--best-root required")
        execute_best(args.best_root,args.output,mu1=args.mu1,
                     grid_n=args.grid_n,vegas_neval=args.neval)
    else:
        print(json.dumps(affinity_protocol(mu1=args.mu1,
                          grid_n=args.grid_n,vegas_neval=args.neval),indent=2))


if __name__=="__main__":
    main()
