"""First actual CLOSED in time χ1↔χ2+massive phi BEST kinetic pilot (no FRW).

Uses independently pinned BEST v1.3.5, full 1<->2 quantum gain/loss,
dynamic mother FD departures, dynamic daughter Bose/Pauli and BEST's
MASSIVE-boson q<qmin zero-mode representation. A finite tiny positive
initial daughter seed 1e-12 avoids an untracked startup seed pass;
the seeded reference number is explicitly included in initial charges.

Gamma_rest is 1 in units T=1, M_chi1=.1, mphi=.001, mchi2=0,
and the bare spin-SUMMED |M|² is adjusted to actual p* as required
by the earlier independent finite-mass source benchmark.

This short fixed-a=1 pilot records two conserved number charges,
sum of PHYSICAL energy including the zero-mode rest energy, IR phi
occupation and zero-mode N0, as well as separate FEASIBILITY/errors.

Do NOT confuse finite qmin Bose population with physical macroscopic
BEC, and do NOT infer FRW frozen primordial 86-eV phase-space f(q).
No external work is done without --execute.
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
from best_massive_daughter_reference import massive_protocol


def dynamic_protocol(mphi=.001,nq=56,vegas_neval=6000,
                     steps=5,dtGamma=1e-5,method="exprb"):
    if steps<1 or dtGamma<=0 or method not in ("exprb","exprb_full","heun"):
        raise ValueError("Positive short evolution steps, duration and supported method required")
    base=massive_protocol(mphi=mphi,grid_n=nq,vegas_neval=vegas_neval)
    return {
        "status":"OPTIONAL CLOSED BEST FINITE-TIME KINETICS; NOT EXECUTED BY DEFAULT",
        "massive_initial_vacuum_protocol":base,
        "daughter_initial_occupation":1e-12,
        "initial_mother":"Fermi-Dirac chi1 mass=.1 at T=1",
        "real_BEST_physics":"time evolving mother+phi+chi2 with Bose/Pauli and detailed balance",
        "static_scale_factor":1.,
        "steps":int(steps),"dtGamma_per_step":float(dtGamma),
        "total_tGamma":float(steps*dtGamma),
        "integrator":method,
        "zero_mode":"BEST finite-mass-boson q<qmin representation; distinguish finite IR cell from true condensate",
        "not_a_physical_86eV_FRW_or_CLASS":True,
        "external_run_executed":False
    }


def execute_best(best_root,output,p):
    folder=Path(best_root).resolve()
    if not (folder/"besthep.py").is_file():
        raise ValueError("BEST besthep.py missing")
    try:
        sha=subprocess.run(["git","-C",str(folder),"rev-parse","HEAD"],
            check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
            text=True,timeout=10).stdout.strip()
    except (OSError,subprocess.SubprocessError) as ex:
        raise ValueError("Cannot verify pinned BEST source") from ex
    if sha!=BEST_PIN:
        raise ValueError("BEST git revision not pinned")
    sys.path.insert(0,str(folder))
    from besthep import BEST,grid_weights
    phys=p["massive_initial_vacuum_protocol"]
    m=phys["mass_chi1"]; mph=phys["mass_phi"]; mchi=phys["mass_chi2"]
    model=BEST(q_min=phys["qmin"],q_max=phys["qmax"],
               n_grid=phys["grid_n"],max_rel_change=.3)
    model.initialize_species("chi1",lambda q:expit(-np.hypot(q,m)),
                             stat="fermion",mass=m,dof=phys["g_chi1"])
    seed=p["daughter_initial_occupation"]
    model.initialize_species("phi",lambda q:seed,
                             stat="boson",mass=mph,dof=phys["g_phi"])
    model.initialize_species("chi2",lambda q:seed,
                             stat="fermion",mass=mchi,dof=phys["g_chi2"])
    bare_m2=phys["spin_summed_fixed_width_matrix_element_squared"]
    model.add_process(
        "dec",["chi1"],["phi","chi2"],
        lambda momenta:np.full(momenta.shape[2],bare_m2),
        neval=phys["vegas_neval_per_grid_point"],nitn=2)
    initial=model.compute_moments()
    t0=float(model.current_time)
    first_charges={
       "chi1_plus_phi":float(initial["chi1"]["n"]+initial["phi"]["n"]),
       "chi1_plus_chi2":float(initial["chi1"]["n"]+initial["chi2"]["n"])
    }
    initial_energy=float(sum(initial[x]["e"] for x in ("chi1","phi","chi2")))
    sequence=[]
    history=model.init_history()
    for step in range(p["steps"]):
        model.evolve_step(p["dtGamma_per_step"],method=p["integrator"],
                          adapt_dt=False)
        moments=model.record(history)
        phi_q=np.asarray(model.r_grids["phi"],float)
        phi_f=np.asarray(model.distributions_1d["phi"],float)
        chi1_f=np.asarray(model.distributions_1d["chi1"],float)
        chi2_f=np.asarray(model.distributions_1d["chi2"],float)
        qmom=grid_weights(phi_q)*phi_q**2*phi_f/(2*math.pi**2)
        numbers={name:float(moments[name]["n"])
                 for name in ("chi1","phi","chi2")}
        energies={name:float(moments[name]["e"])
                  for name in ("chi1","phi","chi2")}
        charge1=numbers["chi1"]+numbers["phi"]
        charge2=numbers["chi1"]+numbers["chi2"]
        total_e=sum(energies.values())
        if min(phi_f.min(),chi1_f.min(),chi2_f.min())<-1e-9:
            raise ArithmeticError("Negative occupation under BEST kinetic time stepping")
        if max(chi1_f.max(),chi2_f.max())>1+1e-8:
            raise ArithmeticError("Pauli bound violated under BEST kinetic time stepping")
        sequence.append({
           "step":step+1,
           "time_after_tGamma":float(model.current_time-t0),
           "number_moments":numbers,
           "physical_energy_moments":energies,
           "total_energy":total_e,
           "charge_chi1_phi_absolute_residual":charge1-first_charges["chi1_plus_phi"],
           "charge_chi1_chi2_absolute_residual":charge2-first_charges["chi1_plus_chi2"],
           "total_physical_energy_absolute_residual":total_e-initial_energy,
           "phi_massive_zero_mode_N0":float(model.zero_mode.get("phi",0.)),
           "phi_smallest_positive_q":float(phi_q[0]),
           "phi_smallest_q_occupation":float(phi_f[0]),
           "phi_maximum_occupation":float(phi_f.max()),
           "chi1_min_occupation":float(chi1_f.min()),
           "chi1_max_occupation":float(chi1_f.max()),
           "chi2_max_occupation":float(chi2_f.max()),
           "physical_phi_q_below_0p01_number_excluding_N0":
               float(qmom[phi_q<.01].sum()),
           "phi_number_including_BEST_zero_mode":numbers["phi"],
           "analytic_initial_vacuum_nphi_rate":
               phys["independent_analytic_phi_number_source"],
           "this_is_not_an_instantaneous_equilibrium_assumption":True
        })
    if model.world_rank!=0:return None
    final=sequence[-1]
    report={
       "external_BEST_executed":True,
       "verified_external_best_git_sha":sha,
       "model":p,
       "initial_number_moments":{name:float(initial[name]["n"]) for name in ("chi1","phi","chi2")},
       "initial_physical_energy_moments":{name:float(initial[name]["e"]) for name in ("chi1","phi","chi2")},
       "initial_two_charges":first_charges,
       "initial_total_energy":initial_energy,
       "time_evolution":sequence,
       "final_time":final["time_after_tGamma"],
       "relative_chi1_phi_charge_error":float(final["charge_chi1_phi_absolute_residual"]/first_charges["chi1_plus_phi"]),
       "relative_chi1_chi2_charge_error":float(final["charge_chi1_chi2_absolute_residual"]/first_charges["chi1_plus_chi2"]),
       "relative_total_physical_energy_error":float((final["total_energy"]-initial_energy)/initial_energy),
       "zero_mode_comes_from_best_massive_phi_true":True,
       "no_Hubble_expansion":True,
       "not_relic_or_class":True,
       "warning":"First short finite-mass closed three-species quantum kinetic time evolution. Source and zero-mode integration still numerical, no FRW/SM plasma, physical 86-eV mass ratio, zero-mode/BEC convergence, CLASS or observations."
    }
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps(report,indent=2))
    return report


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--execute",action="store_true")
    parser.add_argument("--best-root",default="")
    parser.add_argument("--output",default="outputs/best_closed_massive_dynamic.json")
    parser.add_argument("--mphi",type=float,default=.001)
    parser.add_argument("--grid-n",type=int,default=56)
    parser.add_argument("--neval",type=int,default=6000)
    parser.add_argument("--steps",type=int,default=5)
    parser.add_argument("--dt",type=float,default=1e-5)
    parser.add_argument("--method",choices=["exprb","exprb_full","heun"],default="exprb")
    args=parser.parse_args()
    spec=dynamic_protocol(args.mphi,args.grid_n,args.neval,args.steps,
                          args.dt,args.method)
    if args.execute:
        if not args.best_root:raise ValueError("--best-root required")
        execute_best(args.best_root,args.output,spec)
    else:
        print(json.dumps(spec,indent=2))


if __name__=="__main__":
    main()
