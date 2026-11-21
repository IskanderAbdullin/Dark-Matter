"""Reproduce finite-T full-Pauli 2->2 χ-e with arbitrary nonthermal χ.

A deterministic 5D quadrature reference, not independent BEST, not a
physical SM-complete model or cosmological relic result.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import expit
from relativistic_nonlinear_fermion import (
    collision_from_nonthermal_chi, distribution_from_grid)


def audit():
    def f_analytic(eps):
        return lambda q:expit(-np.hypot(q,.1))*(1+eps*np.sin(.7*np.asarray(q)))
    def f_grid(nq,eps,qmax=34.):
        q=np.linspace(0,qmax,nq)
        return distribution_from_grid(q,f_analytic(eps)(q))
    def run(f,n=20,cutoff=30.):
        x=collision_from_nonthermal_chi(
            f,M=.1,Te=1.,G=.001,n_p=n,n_k=n,n_mu=min(16,n),
            n_cm=min(20,n),n_phi=min(16,n),cutoff=cutoff,p_max=cutoff)
        return {key:x[key] for key in (
          "Q_chi_GeV5","Q_e_GeV5","total_entropy_production_GeV4",
          "number_strong_GeV4","energy_strong_GeV5",
          "strong_vs_weak_heat_relative","negative_event_entropy_floor",
          "quadrature_points")}
    cases={
      "analytic_FD_equilibrium":run(f_analytic(0.)),
      "analytic_eps06_n20":run(f_analytic(.06)),
      "analytic_eps06_n24":run(f_analytic(.06),n=24),
      "tabulated_FD_n81_false_signal":run(f_grid(81,0.)),
      "tabulated_FD_n641_false_signal":run(f_grid(641,0.)),
      "tabulated_eps03_n321":run(f_grid(321,.03)),
      "tabulated_eps06_n321":run(f_grid(321,.06)),
      "tabulated_eps12_n321":run(f_grid(321,.12)),
      "tabulated_eps06_n641":run(f_grid(641,.06)),
    }
    assert abs(cases["analytic_FD_equilibrium"]["Q_chi_GeV5"])<1e-18
    a=cases["analytic_eps06_n20"];b=cases["analytic_eps06_n24"]
    assert abs(a["Q_chi_GeV5"]/b["Q_chi_GeV5"]-1)<1e-3
    x=cases["tabulated_FD_n81_false_signal"]
    y=cases["tabulated_FD_n641_false_signal"]
    assert abs(x["Q_chi_GeV5"])>10*abs(y["Q_chi_GeV5"])
    assert x["total_entropy_production_GeV4"]>10*y["total_entropy_production_GeV4"]
    c03=cases["tabulated_eps03_n321"]
    c06=cases["tabulated_eps06_n321"]
    c12=cases["tabulated_eps12_n321"]
    assert abs(c06["Q_chi_GeV5"]/c03["Q_chi_GeV5"]-2)<.02
    assert abs(c12["total_entropy_production_GeV4"]/
               c06["total_entropy_production_GeV4"]-4)<.1
    assert c06["total_entropy_production_GeV4"]>0
    return {
      "status":"deterministic sandbox numerical reference, not cosmological prediction",
      "model":"e+-chi vector contact; Mchi=.1 GeV, Te=1 GeV, G=.001 GeV^-2, mu_e=0",
      "interpolation":"PCHIP 0<=f<=1, hard zero above qmax=34 GeV; artificial heat/entropy expected if under-resolved",
      "benchmarks":cases,
      "does_not_establish":["full cosmological time evolution","SM backreaction",
         "physical gamma_el","electroweak UV completion","86-eV relic",
         "BEC","CLASS-ready"]
    }


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",default="")
    args=parser.parse_args()
    result=audit()
    if args.output:
        path=Path(args.output)
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))
