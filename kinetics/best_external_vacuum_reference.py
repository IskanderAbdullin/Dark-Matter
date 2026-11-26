"""Optional *independent* BEST v1.3.5 quantum 1<->2 reference benchmark.

Source: https://github.com/best-hep/best/tree/c70875be714e7a6bc6d38578c34bb0dc9f5aefa9
BEST is a separate MIT-licensed momentum-resolved Boltzmann solver.
It is NOT vendored or imported by our ordinary tests.

Run only after explicitly installing its optional dependencies:
  git clone https://github.com/best-hep/best.git
  python kinetics/best_external_vacuum_reference.py --protocol
  mpirun -np 2 python kinetics/best_external_vacuum_reference.py \
      --execute --best-root /path/to/best --output outputs/best_vacuum.json

The first stage has mother fFD(E,T=1), EMPTY phi and chi2,
and exact two-body massless daughters. We compare the INTEGRATED
instantaneous rate, not a final nonlinear relic spectrum.
Matrix element convention follows BEST examples/Mdecay.py:
    sum_internal |M|^2 = 16*pi*M*Gamma_rest*g_parent.
Check the independently integrated 3-species rates as well.

When only --protocol is used, NO BEST run is made and NO numerical
agreement is asserted. Third-party Monte Carlo error may be significant;
multiple independent random seeds and q/neval convergence remain required.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import sys
import subprocess

import numpy as np
from scipy.special import expit


BEST_REPOSITORY="best-hep/best"
BEST_PIN="c70875be714e7a6bc6d38578c34bb0dc9f5aefa9"


def benchmark_protocol(mass_over_T=.1,parent_dof=4,rest_width=1.,
                       grid_n=100,qmin=1e-4,qmax=30.,vegas_neval=20000):
    if (mass_over_T<=0 or parent_dof<=0 or rest_width<=0 or grid_n<10
            or qmin<=0 or qmax<=qmin or vegas_neval<1000):
        raise ValueError("Positive physical parameters and workable reference mesh required")
    from onshell_collision_2d import vacuum_reference
    from frw_reservoir_budget import thermal_mother_number
    target=vacuum_reference(M1_T=mass_over_T,g1=parent_dof)
    energy_flux=mass_over_T*rest_width*thermal_mother_number(
        mass_over_T,g1=parent_dof)
    return {
       "status":"PROTOCOL ONLY — no external BEST Monte Carlo was executed",
       "independent_solver":"BEST v1.3.5, MIT-licensed, source best-hep/best",
       "pinned_best_git_commit":BEST_PIN,
       "best_source_url":f"https://github.com/{BEST_REPOSITORY}/tree/{BEST_PIN}",
       "reference_model":"chi1 (fermion, g1=4, mass M/T=.1) -> phi (boson g=1, massless) + chi2 (fermion g=4, massless); initial parent thermal FD T=1 and empty daughters",
       "instantaneous_initial_only":True,
       "T_units":1.,
       "parent_mass_over_T":mass_over_T,
       "rest_width":rest_width,
       "parent_dof":parent_dof,
       "daughter_bose_dof":1,
       "daughter_fermi_dof":4,
       "best_spin_summed_matrix_element_squared":
           float(16*math.pi*mass_over_T*rest_width*parent_dof),
       "grid_n":int(grid_n),"qmin":qmin,"qmax":qmax,
       "vegas_neval_per_grid_point":int(vegas_neval),
       "independently_derived_vacuum_source_n_over_T3_per_time":
           float(target*rest_width),
       "independently_derived_total_daughter_energy_gain_over_T4_per_time":float(energy_flux),
       "analytic_initial_bose_differential_source_q_peak":mass_over_T/2,
       "reference_differential_source":"g1*m*Gamma/(2*pi²)*log(1+exp[-q-m²/(4q)]), vacuum-only",
       "expected_parent_mean_energy_per_decay":float(energy_flux/(target*rest_width)),
       "analytic_energy_identity":"-d(rho_parent)/dt=M1*Gamma_rest*n_parent(T), because Gamma_lab=Gamma_rest*M1/E_parent",
       "continuum_source_normalization":
           "g1*M1*Gamma_rest/(2pi²)*integral_0^inf log(1+exp[-q-M1²/(4q)]) dq",
       "validation_requirements":[
           "T=1 and scale factor a=1; no collision or external bath evolution",
           "constant matrix element convention: summed not averaged over parent helicities",
           "initial daughter occupations numerically zero; no hidden seeding before rate pass",
           "integrate all 3 BEST species instantaneous rates with its own grid_weights",
           "repeat at >=2 grid sizes and >=2 Vegas neval settings; capture Monte Carlo variability",
           "compare parent number loss with BOTH daughter number gains",
           "never call a single snapshot a freeze-out, cold DM or a BEC prediction",
       ],
       "external_run_verified":False
    }


def execute_best(protocol,best_root,output):
    path=Path(best_root).resolve()
    if not (path/"besthep.py").is_file():
        raise ValueError("BEST besthep.py not found at supplied root")
    try:
        actual_sha=subprocess.run(
            ["git","-C",str(path),"rev-parse","HEAD"],check=True,
            stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,
            timeout=10).stdout.strip()
    except (OSError,subprocess.SubprocessError) as error:
        raise ValueError("Independent BEST Git commit provenance unavailable") from error
    if actual_sha!=BEST_PIN:
        raise ValueError(
            "BEST source revision differs from pinned protocol: "+
            actual_sha+" expected "+BEST_PIN)
    sys.path.insert(0,str(path))
    from besthep import BEST, grid_weights  # OPTIONAL external solver
    mass=protocol["parent_mass_over_T"]
    g1=protocol["parent_dof"]
    m2=protocol["best_spin_summed_matrix_element_squared"]
    best=BEST(q_min=protocol["qmin"],q_max=protocol["qmax"],
              n_grid=protocol["grid_n"])
    best.initialize_species(
        "chi1",lambda q:expit(-np.sqrt(np.asarray(q)**2+mass**2)),
        stat="fermion",mass=mass,dof=g1)
    best.initialize_species(
        "phi",lambda q:np.zeros_like(np.asarray(q),dtype=float),
        stat="boson",mass=0.,dof=1)
    best.initialize_species(
        "chi2",lambda q:np.zeros_like(np.asarray(q),dtype=float),
        stat="fermion",mass=0.,dof=4)
    best.add_process(
        "dec",["chi1"],["phi","chi2"],
        lambda momenta:np.full(momenta.shape[2],m2),
        neval=protocol["vegas_neval_per_grid_point"],nitn=2)
    rates=best.compute_rates(["dec"])
    if best.world_rank!=0:
        return None
    rate_integrals={}
    energy_integrals={}
    for species,dof in (("chi1",g1),("phi",1),("chi2",4)):
        q=np.asarray(best.r_grids[species],float)
        collision=np.asarray(rates[species],float)
        rate_integrals[species]=float(
           dof/(2*np.pi**2)*np.dot(grid_weights(q)*q**2,collision))
        on_shell_energy=np.hypot(q,mass) if species=="chi1" else q
        energy_integrals[species]=float(
           dof/(2*np.pi**2)*np.dot(grid_weights(q)*q**2*on_shell_energy,collision))
    from best_vacuum_spectrum import compare_best_phi_grid
    vacuum_shape=compare_best_phi_grid(
        np.asarray(best.r_grids["phi"],float),np.asarray(rates["phi"],float),
        m_over_T=mass,g1=g1,rest_width=protocol["rest_width"],gphi=1)
    produced_phi=rate_integrals["phi"]
    target=protocol["independently_derived_vacuum_source_n_over_T3_per_time"]
    result={
        "external_solver_executed":True,
        "BEST_source_commit_verified":actual_sha,
        "BEST_git_revision_matched_pinned_protocol":True,
        "model_and_protocol":protocol,
        "instantaneous_integrated_number_rates":rate_integrals,
        "instantaneous_integrated_energy_rates":energy_integrals,
        "instantaneous_phi_vacuum_differential_source":vacuum_shape,
        "reference_analytic_parent_energy_loss_magnitude":protocol["independently_derived_total_daughter_energy_gain_over_T4_per_time"],
        "relative_BEST_mother_energy_loss_error":-energy_integrals["chi1"]/protocol["independently_derived_total_daughter_energy_gain_over_T4_per_time"]-1,
        "total_all_species_energy_rate_residual":sum(energy_integrals.values()),
        "reference_analytic_phi_production":target,
        "relative_BEST_phi_rate_error":produced_phi/target-1,
        "chi1_plus_phi_rate_residual":rate_integrals["chi1"]+produced_phi,
        "chi1_plus_chi2_rate_residual":
            rate_integrals["chi1"]+rate_integrals["chi2"],
        "daughter_rates_difference":
            rate_integrals["phi"]-rate_integrals["chi2"],
        "warning":"Single stochastic rate pass, separate q-grid quadrature and possible Monte Carlo fluctuations. Pin external repository yourself; repeat independently before scientific validation. This is not a coupled FRW freeze-out or CLASS input."
    }
    destination=Path(output)
    destination.parent.mkdir(parents=True,exist_ok=True)
    destination.write_text(json.dumps(result,indent=2),encoding="utf-8")
    print(json.dumps(result,indent=2))
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--protocol",action="store_true")
    p.add_argument("--execute",action="store_true")
    p.add_argument("--best-root",default="")
    p.add_argument("--output",default="outputs/best_vacuum.json")
    p.add_argument("--neval",type=int,default=20000)
    p.add_argument("--grid-n",type=int,default=100)
    args=p.parse_args()
    spec=benchmark_protocol(vegas_neval=args.neval,grid_n=args.grid_n)
    if args.execute:
        if not args.best_root:raise ValueError("--best-root required to execute")
        execute_best(spec,args.best_root,args.output)
    else:
        # The safe default is PROTOCOL, not pretending to have run BEST.
        print(json.dumps(spec,indent=2))


if __name__=="__main__":
    main()
