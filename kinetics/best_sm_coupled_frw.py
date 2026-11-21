"""Experimental full BEST 1<->2 massive Bose FRW + explicit SM elastic heat bath.

Physical ingredients that are REALLY coupled here:
  - pinned independent BEST v1.3.5: all chi1<->chi2+massive phi quantum
    gain/loss/Pauli and independently evolving massive bosonic zero mode;
  - prescribed RD a(t)=sqrt(t/t0), comoving q grids, dynamical mother chi1;
  - an elastic, NUMBER-conserving thermal contact with a finite SM
    radiation bath via projected FD(T_SM,mu_chi1) relaxation (BGK);
  - measured dark chi1 heat exchange Q compensated by opposite SM heat
    and feedback into T_SM from a^4 rho_SM = pi²/30*gstar*(aT)^4;
  - pressure-work check for total comoving energy + SM radiation.

Strang-style splitting: elastic half-step at fixed initial a,
BEST n<->m quantum decay full-step, elastic half-step at fixed final a.
Each elastic substep conserves parent number and total dark+SM energy
exactly in its finite quadrature and is Pauli-positivity preserving.
The BEST collision plus the FRW splitting STILL have numerical
time/grid/Monte Carlo errors that MUST be measured, not hidden.

This is a phenomenological elastic gamma_el and ideal constant-gstar
radiation, NOT a lagrangian-derived SM scattering matrix element,
true variable gstarS/gstarRho, Friedmann dark backreaction or
cosmological 86-eV frozen phase-space DM prediction.
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
from best_radiation_frw_pilot import rd_protocol
from sm_radiation_heat_feedback import (
    radiation_elastic_exchange_step,
    comoving_radiation_energy,comoving_radiation_entropy,
)


def sm_coupled_protocol(gamma_el=1e5,gstar=55.,t0=1e-4,
                        steps=8,dtGamma=5e-6,nq=56,neval=6000,
                        method="exprb2"):
    if (not math.isfinite(gamma_el) or gamma_el<0 or
            not math.isfinite(gstar) or gstar<=0):
        raise ValueError("Finite nonnegative elastic rate and positive radiation dof required")
    p=rd_protocol(t0=t0,steps=steps,dtGamma=dtGamma,nq=nq,
                  neval=neval,method=method)
    p["elastic_gamma_el_over_Gamma_rest"]=float(gamma_el)
    p["SM_radiation_gstarS_equals_gstarRho"]=float(gstar)
    p["initial_SM_temperature"]=1.
    p["SM_energy_and_temperature_dynamic"]=True
    p["SM_thermalization_is_number_conserving_BGK_not_microphysical_2to2"]=True
    p["time_splitting"]="elastic half-step, exact BEST quantum 1<->2 full-step, elastic half-step"
    p["FRW_Hubble_prescribed_not_sourced_by_SM_or_dark"]=True
    p["real_dark_matter_86eV_or_CLASS_relic_validated"]=False
    p["external_BEST_executed"]=False
    return p


def execute_best(best_root,output,protocol):
    folder=Path(best_root).resolve()
    if not (folder/"besthep.py").is_file():
        raise ValueError("Independent BEST besthep.py missing")
    try:
        gitsha=subprocess.run(
            ["git","-C",str(folder),"rev-parse","HEAD"],
            check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
            text=True,timeout=10).stdout.strip()
    except (OSError,subprocess.SubprocessError) as e:
        raise ValueError("Pinned external BEST provenance unavailable") from e
    if gitsha!=BEST_PIN:
        raise ValueError("Refusing unpinned independent BEST Git source")
    sys.path.insert(0,str(folder))
    from besthep import BEST,grid_weights,ExtrapolatingInterp

    p=protocol["kinetic_source"]
    M=p["mass_chi1"];mph=p["mass_phi"];mchi=p["mass_chi2"]
    solver=BEST(q_min=p["qmin"],q_max=p["qmax"],
                n_grid=p["grid_n"],max_rel_change=.3)
    solver.set_radiation_dominated(a0=1.,t0=protocol["rd_t0"])
    solver.initialize_species(
        "chi1",lambda q:expit(-np.hypot(q,M)),
        stat="fermion",mass=M,dof=p["g_chi1"])
    solver.initialize_species(
        "phi",lambda q:1e-12*np.exp(-q),
        stat="boson",mass=mph,dof=p["g_phi"])
    solver.initialize_species(
        "chi2",lambda q:1e-12*np.exp(-q),
        stat="fermion",mass=mchi,dof=p["g_chi2"])
    amplitude=p["spin_summed_fixed_width_matrix_element_squared"]
    solver.add_process(
        "dec",["chi1"],["phi","chi2"],
        lambda mom:np.full(mom.shape[2],amplitude),
        neval=p["vegas_neval_per_grid_point"],nitn=2)

    gstar=protocol["SM_radiation_gstarS_equals_gstarRho"]
    ESM=comoving_radiation_energy(1.,1.,gstar)
    gamma_el=protocol["elastic_gamma_el_over_Gamma_rest"]
    elastic_ledger=[]

    def sm_substep(dt,half_label):
        nonlocal ESM
        a=float(solver.scale_factor(solver.current_time))
        T=float((ESM/(math.pi**2/30*gstar*a**4))**.25)
        update=None
        if solver.world_rank==0:
            update=radiation_elastic_exchange_step(
                solver.r_grids["chi1"],solver.distributions_1d["chi1"],
                a=a,mass=M,T_bath=T,g1=p["g_chi1"],
                gamma_el=gamma_el,dt=dt,gstar=gstar)
            solver.distributions_1d["chi1"]=update["mother_after"]
            solver.interpolators["chi1"]=ExtrapolatingInterp(
                solver.r_grids["chi1"],solver.distributions_1d["chi1"],
                mass=M,stat="fermion",a=a)
            ESM=float(update["SM_radiation_comoving_energy_after"])
            entry={
                "label":half_label,
                "time_cosmic":float(solver.current_time),
                "a":a,
                "T_SM_before":T,
                "T_SM_after":float(update["T_SM_after"]),
                "comoving_elastic_dark_energy_delta":
                    float(update["mother_energy_after"]-update["mother_energy_before"]),
                "comoving_elastic_SM_energy_delta":
                    float(update["SM_radiation_comoving_energy_after"]-
                          update["SM_radiation_comoving_energy_before"]),
                "mother_number_residual":float(update["mother_number_change"]),
                "total_entropy_change":float(update["total_entropy_change"]),
                "mu_chi1":float(update["mother_projected_mu"])
            }
        else:
            entry=None
        chi1=solver.world_comm.bcast(
            solver.distributions_1d["chi1"] if solver.world_rank==0 else None,
            root=0)
        interp=solver.world_comm.bcast(
            solver.interpolators["chi1"] if solver.world_rank==0 else None,
            root=0)
        ESM=solver.world_comm.bcast(ESM if solver.world_rank==0 else None,root=0)
        entry=solver.world_comm.bcast(entry,root=0)
        solver.distributions_1d["chi1"]=chi1
        solver.interpolators["chi1"]=interp
        elastic_ledger.append(entry)

    def snapshot():
        a=float(solver.scale_factor(solver.current_time))
        moments=solver.compute_moments()
        darkE=float(sum(moments[name]["e"] for name in ("chi1","phi","chi2")))
        darkP=0.
        for species,deg,mass in (
            ("chi1",p["g_chi1"],M),("phi",p["g_phi"],mph),("chi2",p["g_chi2"],mchi)):
            q=np.asarray(solver.r_grids[species],float)
            f=np.asarray(solver.distributions_1d[species],float)
            E=np.hypot(q,a*mass)
            darkP+=float(deg/(6*np.pi**2)*np.dot(grid_weights(q)*q**4/E,f))
        numbers={name:float(moments[name]["n"]) for name in ("chi1","phi","chi2")}
        totalE=darkE+ESM
        totalP=darkP+ESM/3
        phi_f=np.asarray(solver.distributions_1d["phi"],float)
        return {
            "time":float(solver.current_time),
            "a":a,
            "T_SM":float((ESM/(math.pi**2/30*gstar*a**4))**.25),
            "SM_comoving_energy":float(ESM),
            "SM_comoving_entropy":float(
                comoving_radiation_entropy(a,(ESM/(math.pi**2/30*gstar*a**4))**.25,gstar)),
            "dark_comoving_energy":darkE,
            "dark_comoving_pressure":darkP,
            "total_comoving_energy":totalE,
            "expected_physical_pressure_work_derivative":totalE-3*totalP,
            "number_moments":numbers,
            "charge_chi1_plus_phi":numbers["chi1"]+numbers["phi"],
            "charge_chi1_plus_chi2":numbers["chi1"]+numbers["chi2"],
            "phi_qmin_occupation":float(phi_f[0]),
            "massive_phi_zero_mode_N0":float(solver.zero_mode.get("phi",0.)),
            "chi1_max_occupation":float(
                np.max(solver.distributions_1d["chi1"]))
        }

    initial=snapshot()
    snapshots=[initial]
    pressure_work=0.
    history=solver.init_history()
    for step in range(protocol["steps"]):
        sm_substep(.5*protocol["dtGamma"],"before decay step "+str(step+1))
        solver.evolve_step(protocol["dtGamma"],active_processes=["dec"],
                           method=protocol["integrator"],adapt_dt=False)
        sm_substep(.5*protocol["dtGamma"],"after decay step "+str(step+1))
        solver.record(history)
        new=snapshot()
        previous=snapshots[-1]
        dloga=math.log(new["a"]/previous["a"])
        pressure_work+=.5*dloga*(
            previous["expected_physical_pressure_work_derivative"]+
            new["expected_physical_pressure_work_derivative"])
        new["accumulated_pressure_work"]=pressure_work
        new["total_energy_change_from_initial"]=(
            new["total_comoving_energy"]-initial["total_comoving_energy"])
        new["physical_energy_closure_residual"]=(
            new["total_energy_change_from_initial"]-pressure_work)
        new["charge_chi1_plus_phi_residual"]=(
            new["charge_chi1_plus_phi"]-initial["charge_chi1_plus_phi"])
        new["charge_chi1_plus_chi2_residual"]=(
            new["charge_chi1_plus_chi2"]-initial["charge_chi1_plus_chi2"])
        snapshots.append(new)

    if solver.world_rank!=0:return None
    last=snapshots[-1]
    totalE0=initial["total_comoving_energy"]
    charge0=initial["charge_chi1_plus_phi"]
    out={
        "independent_BEST_and_finite_SM_bath_ACTUALLY_coupled":True,
        "verified_external_BEST_git_sha":gitsha,
        "protocol":protocol,
        "initial_snapshot":initial,
        "time_series":snapshots,
        "elastic_heat_exchange_substeps":elastic_ledger,
        "final_T_SM":last["T_SM"],
        "net_SM_comoving_energy_change":
            last["SM_comoving_energy"]-initial["SM_comoving_energy"],
        "net_dark_comoving_energy_change":
            last["dark_comoving_energy"]-initial["dark_comoving_energy"],
        "final_relative_total_energy_minus_pressure_work":
            last["physical_energy_closure_residual"]/totalE0,
        "final_relative_charge_chi1_phi_error":
            last["charge_chi1_plus_phi_residual"]/charge0,
        "final_relative_charge_chi1_chi2_error":
            last["charge_chi1_plus_chi2_residual"]/charge0,
        "time_step_convergence_confirmed":False,
        "SM_gstarS_temperature_dependence_included":False,
        "actual_SM_2to2_matrix_element_included":False,
        "Friedmann_dark_backreaction_included":False,
        "physical_86eV_DM_or_CLASS_ready":False,
        "warning":"True coupled BEST+finite ideal-radiation heat Q, but gamma_el is phenomenological and only elastic chi1-number-conserving contact. One short prescribed-RD run does not validate primordial eV dark matter."
    }
    file=Path(output);file.parent.mkdir(parents=True,exist_ok=True)
    file.write_text(json.dumps(out,indent=2),encoding="utf-8")
    print(json.dumps(out,indent=2))
    return out


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--execute",action="store_true")
    parser.add_argument("--best-root",default="")
    parser.add_argument("--output",default="outputs/best_sm_coupled_frw.json")
    parser.add_argument("--gamma-el",type=float,default=1e5)
    parser.add_argument("--gstar",type=float,default=55.)
    parser.add_argument("--steps",type=int,default=8)
    parser.add_argument("--dt",type=float,default=5e-6)
    parser.add_argument("--grid-n",type=int,default=56)
    parser.add_argument("--neval",type=int,default=6000)
    parser.add_argument("--method",default="exprb2",choices=["exprb","exprb2"])
    opts=parser.parse_args()
    p=sm_coupled_protocol(opts.gamma_el,opts.gstar,
                          steps=opts.steps,dtGamma=opts.dt,
                          nq=opts.grid_n,neval=opts.neval,method=opts.method)
    if opts.execute:
        if not opts.best_root:raise ValueError("--best-root required")
        execute_best(opts.best_root,opts.output,p)
    else:
        print(json.dumps(p,indent=2))


if __name__=="__main__":
    main()
