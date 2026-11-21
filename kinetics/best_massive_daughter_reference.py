"""Optional independent BEST vacuum decay of a MASSIVE Bose daughter.

At T=1, a thermal FD mother M1=.1 and initial empty phi(mphi>0),
chi2(mchi>=0). Gamma_rest is held FIXED while varying daughter masses.
Therefore the necessary bare spin-SUMMED matrix element is

  sum |M|² = g1 * 8*pi*M1²*Gamma_rest/p_star

where p_star is the physical rest-frame daughter momentum.
The massless limit is g1*16*pi*M1*Gamma_rest, consistent with the
already independently validated BEST vacuum source.

References: massive_daughter_decay_source.py (analytic spectrum, IR)
and the pinned independent BEST c70875be714e7a6bc6d38578c34bb0dc9f5aefa9.

Numerical side: BEST's finite massive boson can have a zero-mode
representation; ordinary phi momentum rates are compared only on its
resolved positive-q grid, NOT used to infer n0 or an actual condensate.
No external computation occurs without --execute.
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
from massive_daughter_decay_source import (
    rest_kinematics, integrated_source_and_energy,
    massive_daughter_vacuum_spectrum,
    zero_momentum_bose_occupation_source,
)


def massive_protocol(mass=.1,mphi=.001,mchi=0.,g1=4,rest_width=1.,
                     qmin=1e-4,qmax=30.,grid_n=56,vegas_neval=12000):
    k=rest_kinematics(mass,mphi,mchi)
    if (g1<=0 or rest_width<=0 or qmin<=0 or qmax<=qmin or
            grid_n<10 or vegas_neval<1000):
        raise ValueError("Positive physical parameters and workable q grid required")
    analytic=integrated_source_and_energy(mass,mphi,mchi,g1,rest_width)
    matrix_squared=g1*8*math.pi*mass**2*rest_width/k["rest_daughter_momentum"]
    return {
        "status":"ANALYTIC MASSIVE DAUGHTER PROTOCOL; EXTERNAL BEST NOT EXECUTED",
        "verified_best_git_commit_expected":BEST_PIN,
        "mass_chi1":mass,"mass_phi":mphi,"mass_chi2":mchi,
        "g_chi1":g1,"g_phi":1,"g_chi2":4,
        "Gamma_rest":rest_width,"T":1.,
        "qmin":qmin,"qmax":qmax,"grid_n":grid_n,
        "vegas_neval_per_grid_point":vegas_neval,
        "daughter_rest_momentum":k["rest_daughter_momentum"],
        "spin_summed_fixed_width_matrix_element_squared":matrix_squared,
        "independent_analytic_phi_number_source":analytic["fixed_gamma_exact_number_rate"],
        "independent_analytic_parent_energy_loss":analytic["fixed_gamma_parent_physical_energy_loss"],
        "independent_analytic_phi_energy_rate":analytic["rest_frame_exact_phi_physical_energy_rate"],
        "independent_analytic_chi2_energy_rate":
            analytic["fixed_gamma_parent_physical_energy_loss"]*
            analytic["exact_chi_energy_fraction"],
        "independent_exact_q0_phi_occupation_rate":
            analytic["finite_mass_q0_Bose_occupation_source"],
        "initial_daughters_empty":True,
        "assumed_fixed_rest_width_not_fixed_matrix_element":True,
        "not_physical_86eV_BEC_or_CLASS":True
    }


def execute_best(protocol,best_root,output):
    path=Path(best_root).resolve()
    if not (path/"besthep.py").is_file():
        raise ValueError("Independent BEST besthep.py missing")
    try:
        sha=subprocess.run(
            ["git","-C",str(path),"rev-parse","HEAD"],
            check=True,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
            timeout=10).stdout.strip()
    except (OSError,subprocess.SubprocessError) as error:
        raise ValueError("Independent BEST Git provenance unavailable") from error
    if sha!=BEST_PIN:
        raise ValueError("BEST Git SHA mismatch; refusing unpinned science")
    sys.path.insert(0,str(path))
    from besthep import BEST,grid_weights
    m=protocol["mass_chi1"]
    mph=protocol["mass_phi"]
    mch=protocol["mass_chi2"]
    model=BEST(q_min=protocol["qmin"],q_max=protocol["qmax"],
               n_grid=protocol["grid_n"])
    model.initialize_species(
        "chi1",lambda q:expit(-np.hypot(q,m)),
        stat="fermion",mass=m,dof=protocol["g_chi1"])
    model.initialize_species(
        "phi",lambda q:np.zeros_like(np.asarray(q),dtype=float),
        stat="boson",mass=mph,dof=protocol["g_phi"])
    model.initialize_species(
        "chi2",lambda q:np.zeros_like(np.asarray(q),dtype=float),
        stat="fermion",mass=mch,dof=protocol["g_chi2"])
    matrix_element=protocol["spin_summed_fixed_width_matrix_element_squared"]
    model.add_process(
        "dec",["chi1"],["phi","chi2"],
        lambda mom:np.full(mom.shape[2],matrix_element),
        neval=protocol["vegas_neval_per_grid_point"],nitn=2)
    rate=model.compute_rates(["dec"])
    if model.world_rank!=0:
        return None
    sources={}
    for species,deg,mass in (("chi1",protocol["g_chi1"],m),
                             ("phi",protocol["g_phi"],mph),
                             ("chi2",protocol["g_chi2"],mch)):
        q=np.asarray(model.r_grids[species],dtype=float)
        C=np.asarray(rate[species],dtype=float)
        vol=deg/(2*np.pi**2)*grid_weights(q)*q*q
        sources[species]={
            "number_rate":float(np.dot(vol,C)),
            "physical_energy_rate":float(np.dot(vol*np.hypot(q,mass),C)),
            "min_occupation_source":float(np.min(C)),
            "max_occupation_source":float(np.max(C))
        }
    q=np.asarray(model.r_grids["phi"],float)
    phi_rate=np.asarray(rate["phi"],float)
    observed=protocol["g_phi"]/(2*np.pi**2)*q*q*phi_rate
    predicted=massive_daughter_vacuum_spectrum(
        q,M=m,mphi=mph,mchi=mch,g1=protocol["g_chi1"],
        rest_width=protocol["Gamma_rest"])
    peak=float(np.max(predicted))
    physical=predicted>peak*.001
    if peak<=0 or np.count_nonzero(physical)<3:
        raise ValueError("Independent positive-q BEST mesh does not cover massive source")
    rms=float(np.sqrt(np.mean(((observed[physical]-predicted[physical])/peak)**2)))
    reference_number=protocol["independent_analytic_phi_number_source"]
    parent_energy=protocol["independent_analytic_parent_energy_loss"]
    result={
      "external_BEST_actually_executed":True,
      "actual_git_sha":sha,
      "protocol":protocol,
      "positive_q_resolved_species_rates":sources,
      "real_phi_source_q":q.tolist(),
      "real_phi_df_dt":phi_rate.tolist(),
      "analytic_phi_dR_dq":predicted.tolist(),
      "best_phi_dR_dq":observed.tolist(),
      "relative_phi_number_rate_error":
          sources["phi"]["number_rate"]/reference_number-1,
      "relative_mother_energy_loss_error":
          -sources["chi1"]["physical_energy_rate"]/parent_energy-1,
      "relative_phi_energy_rate_error":
          sources["phi"]["physical_energy_rate"]/
          protocol["independent_analytic_phi_energy_rate"]-1,
      "total_three_species_physical_energy_residual":
          float(sum(item["physical_energy_rate"] for item in sources.values())),
      "relative_three_species_number_residual":
          (sources["chi1"]["number_rate"]+sources["phi"]["number_rate"])/reference_number,
      "positive_q_shape_peak_normalized_rms_error":rms,
      "smallest_q":float(q[0]),
      "best_df_phi_dt_at_smallest_positive_q":float(phi_rate[0]),
      "analytic_df_phi_dt_at_smallest_positive_q":
          float(2*np.pi**2*predicted[0]/(protocol["g_phi"]*q[0]**2)),
      "analytic_true_q0_occupation_source":
          protocol["independent_exact_q0_phi_occupation_rate"],
      "zero_mode_initial_occupation_assumed":0.,
      "condensate_production_or_dynamics_validated":False,
      "CLASS_ready":False,
      "warning":"Fixed rest Gamma, finite-mass two-body vacuum ONLY. The BEST qmin mode is not a homogeneous condensate. No closed dynamical Bose/FRW or physical 86-eV relic inferred."
    }
    dest=Path(output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,indent=2),encoding="utf-8")
    print(json.dumps(result,indent=2))
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--execute",action="store_true")
    p.add_argument("--best-root",default="")
    p.add_argument("--output",default="outputs/best_massive_phi.json")
    p.add_argument("--mphi",type=float,default=.001)
    p.add_argument("--mchi",type=float,default=0.)
    p.add_argument("--grid-n",type=int,default=56)
    p.add_argument("--neval",type=int,default=12000)
    opts=p.parse_args()
    model=massive_protocol(mphi=opts.mphi,mchi=opts.mchi,
                           grid_n=opts.grid_n,vegas_neval=opts.neval)
    if opts.execute:
        if not opts.best_root:raise ValueError("--best-root required")
        execute_best(model,opts.best_root,opts.output)
    else:
        print(json.dumps(model,indent=2))


if __name__=="__main__":
    main()
