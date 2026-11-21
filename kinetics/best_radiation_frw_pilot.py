"""First independent full quantum 3-species BEST evolution in prescribed RD FRW.

Unlike the previous fixed-a toy, the external BEST solver now uses
  a(t) = sqrt(t/t0) for t>=t0, a(t0)=1
with massive chi1, massive Bose phi and massless Fermi chi2 on
FIXED COMOVING momentum q grids. All 1<->2 gain/loss channels and
the massive phi zero mode evolve dynamically, without resetting the
mother to instantaneous thermal FD.

Independent physical FRW audit:
 BEST comoving E = sum[g/(2pi²) integral dq q² f sqrt(q²+(a m)²)
                      + a*m*N0].
 Define Pcom = sum[g/(6pi²) integral dq q⁴ f/Ecom].
 For a CLOSED collision network in prescribed FRW:
 dEcom/dln a = Ecom - 3 Pcom.
 The pressure work is numerically integrated using the actual
 distributions on both ends of each time step (trapezoid in ln a).
 Both comoving charges Nchi1+Nphi, Nchi1+Nchi2 must stay constant.

This tests a RADIATION BACKGROUND prescribed externally: SM entropy
degrees fixed, dark-sector backreaction neglected, masses and Gamma
constants, no physical Γ(T) or fully realistic cosmological history.
The initial state is an arbitrary thermal FD mother at T=1 with
daughter seeds 1e-12 exp(-q); not a primordial DM relic prediction.
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
from best_massive_daughter_reference import massive_protocol
from best_external_vacuum_reference import BEST_PIN


def rd_protocol(t0=.0001,steps=4,dtGamma=1e-5,nq=56,
                neval=6000,phi_mass=.001,method="exprb"):
    if (t0<=0 or steps<1 or dtGamma<=0 or
            method not in ("exprb","exprb2","exprb_full","heun")):
        raise ValueError("Positive radiation cosmic time, step count and supported method required")
    kinetic=massive_protocol(mphi=phi_mass,grid_n=nq,
                            vegas_neval=neval)
    aend=math.sqrt(1+steps*dtGamma/t0)
    return {
        "status":"CONDITIONAL PRESCRIBED-RD BEST 3-SPECIES FRW PROTOCOL; NOT EXECUTED BY DEFAULT",
        "kinetic_source":kinetic,
        "rd_t0":float(t0),
        "rd_initial_scale_factor":1.,
        "rd_final_expected_a":aend,
        "rd_initial_H_over_Gamma":1/(2*t0),
        "fixed_SM_entropy_degrees":True,
        "constant_physical_masses_and_Gamma":True,
        "steps":int(steps),"dtGamma":float(dtGamma),
        "end_tGamma":float(steps*dtGamma),
        "integrator":method,
        "daughter_initial_seed":"1e-12 exp(-q)",
        "dynamic_mother_is_NOT_instantaneously_rethermalized":True,
        "full_SM_energy_transfer_or_entropy_included":False,
        "self_consistent_Hubble_from_dark_sector":False,
        "primordial_86eV_relic_or_CLASS_ready":False
    }


def execute_best(best_root,output,p):
    folder=Path(best_root).resolve()
    if not (folder/"besthep.py").is_file():
        raise ValueError("Independent BEST besthep.py missing")
    try:
        sha=subprocess.run(["git","-C",str(folder),"rev-parse","HEAD"],
            check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
            text=True,timeout=10).stdout.strip()
    except (OSError,subprocess.SubprocessError) as ex:
        raise ValueError("Cannot validate BEST Git SHA") from ex
    if sha!=BEST_PIN:
        raise ValueError("BEST Git SHA differs from pinned external source")
    sys.path.insert(0,str(folder))
    from besthep import BEST,grid_weights
    kin=p["kinetic_source"]
    M=kin["mass_chi1"];mphi=kin["mass_phi"];mchi=kin["mass_chi2"]
    g1=kin["g_chi1"];gphi=kin["g_phi"];gchi=kin["g_chi2"]
    solver=BEST(q_min=kin["qmin"],q_max=kin["qmax"],
                n_grid=kin["grid_n"],max_rel_change=.3)
    solver.set_radiation_dominated(a0=1.,t0=p["rd_t0"])
    solver.initialize_species(
        "chi1",lambda q:expit(-np.hypot(q,M)),
        stat="fermion",mass=M,dof=g1)
    solver.initialize_species(
        "phi",lambda q:1e-12*np.exp(-q),
        stat="boson",mass=mphi,dof=gphi)
    solver.initialize_species(
        "chi2",lambda q:1e-12*np.exp(-q),
        stat="fermion",mass=mchi,dof=gchi)
    amp=kin["spin_summed_fixed_width_matrix_element_squared"]
    solver.add_process(
        "dec",["chi1"],["phi","chi2"],
        lambda mom:np.full(mom.shape[2],amp),
        neval=kin["vegas_neval_per_grid_point"],nitn=2)

    def snapshot():
        moments=solver.compute_moments()
        a=solver.scale_factor(solver.current_time)
        pressure={}
        charges={}
        for name,deg,mass in (("chi1",g1,M),("phi",gphi,mphi),("chi2",gchi,mchi)):
            q=np.asarray(solver.r_grids[name],dtype=float)
            f=np.asarray(solver.distributions_1d[name],dtype=float)
            E=np.hypot(q,a*mass)
            pressure[name]=float(deg/(6*np.pi**2)*np.dot(
                grid_weights(q)*q**4/E,f))
        total_e=float(sum(moments[n]["e"] for n in ("chi1","phi","chi2")))
        total_p=sum(pressure.values())
        charges["chi1_plus_phi"]=float(moments["chi1"]["n"]+moments["phi"]["n"])
        charges["chi1_plus_chi2"]=float(moments["chi1"]["n"]+moments["chi2"]["n"])
        fq=np.asarray(solver.distributions_1d["phi"],dtype=float)
        return {
            "time_cosmic":float(solver.current_time),
            "time_since_initial_Gamma":float(solver.current_time-p["rd_t0"]),
            "a":float(a),"x_M1_over_T":float(M*a),
            "H_over_Gamma":float(1/(2*solver.current_time)),
            "species_comoving_numbers":{k:float(moments[k]["n"]) for k in ("chi1","phi","chi2")},
            "species_comoving_energies":{k:float(moments[k]["e"]) for k in ("chi1","phi","chi2")},
            "species_comoving_pressures":pressure,
            "total_comoving_energy":total_e,
            "total_comoving_pressure":total_p,
            "expected_dEcom_dln_a_for_closed_source":total_e-3*total_p,
            "two_number_charges":charges,
            "phi_massive_zero_mode_N0":float(solver.zero_mode.get("phi",0.)),
            "phi_occupation_at_smallest_positive_q":float(fq[0]),
            "phi_maximum_occupation":float(fq.max()),
            "not_instantaneous_FD_bath":True
        }

    first=snapshot()
    history=solver.init_history()
    timeseries=[first]
    pressure_work=0.
    for step in range(p["steps"]):
        solver.evolve_step(p["dtGamma"],method=p["integrator"],adapt_dt=False)
        solver.record(history)
        now=snapshot()
        previous=timeseries[-1]
        dx=math.log(now["a"]/previous["a"])
        pressure_work+=.5*dx*(
            previous["expected_dEcom_dln_a_for_closed_source"]+
            now["expected_dEcom_dln_a_for_closed_source"])
        now["independent_accumulated_FRW_pressure_work"]=pressure_work
        now["comoving_energy_change_from_initial"] = (
            now["total_comoving_energy"]-first["total_comoving_energy"])
        now["physical_energy_closure_residual"]= (
            now["comoving_energy_change_from_initial"]-pressure_work)
        now["chi1_plus_phi_charge_residual"]= (
            now["two_number_charges"]["chi1_plus_phi"]-
            first["two_number_charges"]["chi1_plus_phi"])
        now["chi1_plus_chi2_charge_residual"]= (
            now["two_number_charges"]["chi1_plus_chi2"]-
            first["two_number_charges"]["chi1_plus_chi2"])
        timeseries.append(now)
    if solver.world_rank!=0:return None
    last=timeseries[-1]
    norm=first["total_comoving_energy"]
    charge_scale=first["two_number_charges"]["chi1_plus_phi"]
    data={
        "real_independent_BEST_FRW_executed":True,
        "verified_best_git_sha":sha,
        "model":p,
        "initial":first,"time_series":timeseries,
        "accumulated_expected_radiation_FRW_pressure_work":
            pressure_work,
        "observed_total_comoving_energy_change":
            last["comoving_energy_change_from_initial"],
        "independent_physical_energy_closure_residual":
            last["physical_energy_closure_residual"],
        "relative_closed_FRW_energy_error":
            last["physical_energy_closure_residual"]/norm,
        "relative_charge_chi1_phi_error":
            last["chi1_plus_phi_charge_residual"]/charge_scale,
        "relative_charge_chi1_chi2_error":
            last["chi1_plus_chi2_charge_residual"]/charge_scale,
        "massive_zero_mode_actually_evolved":True,
        "SM_bath_energy_Q_simulated":False,
        "Hubble_dark_sector_backreaction_simulated":False,
        "zero_mode_macroscopic_BEC_proven":False,
        "frozen_dark_matter_PSD_or_CLASS_ready":False,
        "warning":"Prescribed RD background, constant gstarS and Γ, short time only, possible MC/time discretization energy error. Zero mode represents BEST q<qmin basis and is not BEC evidence."
    }
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(data,indent=2),encoding="utf-8")
    print(json.dumps(data,indent=2))
    return data


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--execute",action="store_true")
    ap.add_argument("--best-root",default="")
    ap.add_argument("--output",default="outputs/best_prescribed_rd_frw.json")
    ap.add_argument("--t0",type=float,default=.0001)
    ap.add_argument("--steps",type=int,default=4)
    ap.add_argument("--dt",type=float,default=1e-5)
    ap.add_argument("--grid-n",type=int,default=56)
    ap.add_argument("--neval",type=int,default=6000)
    ap.add_argument("--method",choices=["exprb","exprb2","exprb_full","heun"],default="exprb")
    o=ap.parse_args()
    p=rd_protocol(o.t0,o.steps,o.dt,o.grid_n,o.neval,method=o.method)
    if o.execute:
        if not o.best_root:raise ValueError("--best-root required")
        execute_best(o.best_root,o.output,p)
    else:
        print(json.dumps(p,indent=2))


if __name__=="__main__":
    main()
