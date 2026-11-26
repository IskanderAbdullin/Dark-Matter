"""Independent 1<->2 quantum detailed-balance benchmark for BEST v1.3.5.

At the same T, parent FD and daughter FD/BE occupations with
  mu_chi1 = mu_chi2 + mu_phi
obey EXACT microscopic forward/reverse balance on every on-shell event:
  f1*(1-f2)*(1+bphi) = (1-f1)*f2*bphi
provided E1=E2+Ephi.

Default mu_phi=-.3 (finite boson IR), mu_chi2=.1, mu_chi1=-.2.
The optional external BEST run computes all three mode-resolved net,
physical gain/loss and integrated number+energy rates. The default
operation merely checks the analytic event identity and NEVER claims
BEST was executed.

The source is independently pinned to BEST c70875be714e7a6bc6d38578c34bb0dc9f5aefa9
and rejects any mismatched Git checkout. No external packages are
needed by the ordinary in-repository regression tests.
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

from best_external_vacuum_reference import BEST_PIN, benchmark_protocol


def fermi(E,mu=0.,temperature=1.):
    if temperature<=0:
        raise ValueError("positive T")
    return expit(-(np.asarray(E)-mu)/temperature)


def bose(E,mu=-.3,temperature=1.):
    if temperature<=0 or np.any(np.asarray(E)-mu<=0):
        raise ValueError("Bose state requires E>mu and positive T")
    return 1/np.expm1((np.asarray(E)-mu)/temperature)


def event_balance(qphi,qchi,m=.1,mu_phi=-.3,mu_chi=.1,
                  mu_parent=None):
    if m<=0 or qphi<=0 or qchi<=0 or 4*qphi*qchi<m*m:
        raise ValueError("Positive physical on-shell daughter energies required")
    parent_mu=mu_phi+mu_chi if mu_parent is None else mu_parent
    f1=float(fermi(qphi+qchi,parent_mu))
    f2=float(fermi(qchi,mu_chi))
    b=float(bose(qphi,mu_phi))
    A=f1*(1-f2)*(1+b)
    B=(1-f1)*f2*b
    return {
        "qphi":float(qphi),"qchi":float(qchi),
        "on_shell_parent_energy":float(qphi+qchi),
        "mu_parent":parent_mu,"mu_chi2":mu_chi,"mu_phi":mu_phi,
        "forward_quantum_factor":A,
        "inverse_quantum_factor":B,
        "absolute_gain_loss_residual":A-B,
        "relative_balance_error":(A-B)/max(A,B,1e-300),
        "chemical_equilibrium":abs(parent_mu-mu_phi-mu_chi)<1e-12
    }


def protocol():
    vacuum=benchmark_protocol()
    mu_phi=-.3;mu_chi=.1;mu_parent=mu_phi+mu_chi
    grid=np.geomspace(.005,2.0,40)
    values=[]
    for x in grid:
        y=max(.1*.1/(4*x)*(1+1e-12),.005)
        values.append(event_balance(float(x),float(y),
                                    mu_phi=mu_phi,mu_chi=mu_chi))
    max_error=max(abs(row["relative_balance_error"]) for row in values)
    if max_error>1e-11:
        raise AssertionError("Analytic Bose-Fermi detailed balance failed")
    return {
        "status":"ANALYTIC EQUILIBRIUM ONLY; EXTERNAL BEST NOT EXECUTED",
        "best_pin":BEST_PIN,
        "mass_over_T":.1,
        "Gamma_rest":1.,
        "g1":4,"gphi":1,"gchi2":4,
        "mu_parent":mu_parent,"mu_phi":mu_phi,"mu_chi2":mu_chi,
        "best_spin_summed_amplitude":vacuum["best_spin_summed_matrix_element_squared"],
        "max_pointwise_analytic_gain_loss_error":max_error,
        "quantum_statistics":"parent fermion, daughter phi boson, daughter chi2 fermion",
        "required_identity":"mu_chi1=mu_phi+mu_chi2; E_chi1=E_phi+E_chi2",
        "number_of_independent_analytic_events":len(values),
        "external_best_executed":False,
        "not_cosmological_freezeout":True
    }


def execute_best(best_root,output,neval=20000,grid_n=100):
    spec=protocol()
    path=Path(best_root).resolve()
    if not (path/"besthep.py").is_file():
        raise ValueError("BEST besthep.py not found")
    try:
        sha=subprocess.run(["git","-C",str(path),"rev-parse","HEAD"],
            text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
            check=True,timeout=10).stdout.strip()
    except (OSError,subprocess.SubprocessError) as error:
        raise ValueError("Cannot verify BEST Git provenance") from error
    if sha!=BEST_PIN:
        raise ValueError("BEST Git SHA mismatch: "+sha)
    sys.path.insert(0,str(path))
    from besthep import BEST,grid_weights
    best=BEST(q_min=1e-4,q_max=30.,n_grid=grid_n)
    mp=spec["mass_over_T"]
    best.initialize_species(
        "chi1",lambda q:fermi(np.hypot(q,mp),spec["mu_parent"]),
        stat="fermion",mass=mp,dof=4)
    best.initialize_species(
        "phi",lambda q:bose(np.asarray(q),spec["mu_phi"]),
        stat="boson",mass=0.,dof=1)
    best.initialize_species(
        "chi2",lambda q:fermi(np.asarray(q),spec["mu_chi2"]),
        stat="fermion",mass=0.,dof=4)
    mesq=spec["best_spin_summed_amplitude"]
    best.add_process(
        "dec",["chi1"],["phi","chi2"],
        lambda momenta:np.full(momenta.shape[2],mesq),
        neval=neval,nitn=2)
    raw=best.compute_rates(["dec"])
    if best.world_rank!=0:
        return None
    by_species={}
    for species,dof,mass in (("chi1",4,mp),("phi",1,0.),("chi2",4,0.)):
        q=np.asarray(best.r_grids[species],float)
        weight=grid_weights(q)*q*q*dof/(2*np.pi**2)
        net=np.asarray(raw[species],float)
        gain=np.asarray(best.gain_rates[species],float)
        loss=np.asarray(best.loss_rates[species],float)
        gross=np.dot(weight,gain+loss)
        net_num=np.dot(weight,net)
        energy=np.hypot(q,mass)
        net_energy=np.dot(weight*energy,net)
        local_scale=np.maximum(gain,loss)
        visible=local_scale>max(float(np.max(local_scale))*.001,1e-300)
        by_species[species]={
            "number_rate":float(net_num),
            "number_relative_to_gross_gain_loss":float(
                abs(net_num)/max(float(gross),1e-300)),
            "energy_rate":float(net_energy),
            "maximum_local_net_rate_relative_to_gross_on_support":
                float(np.max(abs(net[visible])/local_scale[visible]))
                if np.any(visible) else None,
            "mode_count":int(len(q))
        }
    result={
        "external_solver_executed":True,
        "best_actual_verified_git_SHA":sha,
        "protocol":spec,
        "q_grid_n":grid_n,"vegas_neval":neval,
        "all_three_species_equilibrium":by_species,
        "net_number_charge_residuals":{
            "parent_plus_phi":
               by_species["chi1"]["number_rate"]+by_species["phi"]["number_rate"],
            "parent_plus_chi2":
               by_species["chi1"]["number_rate"]+by_species["chi2"]["number_rate"]},
        "total_net_energy_rate":sum(x["energy_rate"] for x in by_species.values()),
        "warning":"An actual stochastic BEST equilibrium pass; verify multiple seeds, neval and grids. Not the cold nonlinear Bose burst, condensate or physical FRW relic."
    }
    dest=Path(output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,indent=2),encoding="utf-8")
    print(json.dumps(result,indent=2))
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--execute",action="store_true")
    p.add_argument("--best-root",default="")
    p.add_argument("--grid-n",type=int,default=100)
    p.add_argument("--neval",type=int,default=20000)
    p.add_argument("--output",default="outputs/best_thermal_balance.json")
    args=p.parse_args()
    if args.execute:
        if not args.best_root:raise ValueError("--best-root required")
        execute_best(args.best_root,args.output,args.neval,args.grid_n)
    else:
        print(json.dumps(protocol(),indent=2))


if __name__=="__main__":
    main()
