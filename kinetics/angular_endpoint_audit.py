"""Angular endpoint resolution audit for massless daughters of M1/T=0.1.

For isotropic two-body decay, q=(sqrt(p1^2+M1^2)+p1*mu)/2.
Gauss-Legendre angular quadrature has NO mu=-1 node. This creates an
artificial infrared support floor above the physical continuum endpoint.
The same issue affects the original event kernel and the lattice surrogate.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path

from scipy.special import roots_legendre
from one_to_two_fixed_bath import FlatThermalDecay
from lattice_detailed_balance import EntropyBalancedLattice


def analytic_angular_floor(nmu, M1_T=0.1, p_max=26.):
    if nmu<2 or M1_T<=0 or p_max<=0:
        raise ValueError("Invalid quadrature or physical parameters")
    mu=float(roots_legendre(nmu)[0][0])
    # mu<0. With fixed angular node, minimize q(p,mu) in parent momentum.
    p_stationary=M1_T*(-mu)/math.sqrt(1-mu*mu)
    p_star=min(p_max,max(0.,p_stationary))
    floor=(math.hypot(p_star,M1_T)+p_star*mu)/2
    true_floor=(math.hypot(p_max,M1_T)-p_max)/2
    return {
        "nmu":nmu,
        "mu_most_backward":mu,
        "parent_p_optimal_at_mu":p_stationary,
        "parent_p_used_after_pmax_bound":p_star,
        "analytic_angular_support_floor_q":floor,
        "true_continuum_angle_floor_for_p_le_26":true_floor,
        "angular_floor_over_continuum_floor":floor/true_floor,
    }


def quadrature_support(nmu,np1=240,nq=150,M1_T=.1):
    k=FlatThermalDecay(nq=nq,np1=np1,nmu=nmu,M1_T=M1_T)
    angular=analytic_angular_floor(nmu,M1_T=M1_T)
    observed=float(min(k.qphi.min(),k.qchi.min()))
    angular.update({
        "np1":np1, "nq":nq,
        "lowest_actual_quadrature_daughter_q":observed,
        "relative_excess_from_finite_parent_quadrature":observed/angular["analytic_angular_support_floor_q"]-1,
        "nominal_q_grid_min":float(k.edges[0]),
        "angular_support_floor_over_nominal_qmin":observed/k.edges[0],
        "qmin_below_all_generated_daughter_events":bool(k.edges[0]<observed)
    })
    return angular


def angular_spectrum(nmu,np1=240,nq=220):
    m=EntropyBalancedLattice(nq=nq,np1=np1,nmu=nmu,qmin=1e-5)
    evolved=m.evolve()
    return {**quadrature_support(nmu,np1,nq),
            "surrogate_evolved":evolved,
            "surrogate_vacuum_rate_relative_shift":m.relative_vacuum_shift}


def run(quick=False):
    angular_nodes=(48,72,120) if quick else (48,72,96,120,160,200,300,400)
    cases=[angular_spectrum(nmu) for nmu in angular_nodes]
    for row in cases:
        assert row["qmin_below_all_generated_daughter_events"],row
        assert row["relative_excess_from_finite_parent_quadrature"]>=-1e-10,row
    return {
        "scope":"fixed mother, physical 2-body angular quadrature; lattice-SURROGATE finite-time spectrum",
        "physics":"M1/T=0.1, massless daughters, mother p1 integration truncated at p/T=26",
        "angular_quadrature":"Gauss-Legendre on mu in [-1,1], omits endpoints",
        "not_a_cosmological_spectrum":True,
        "results":cases,
        "main_conclusion":"Changing nominal qmin below the smallest angular-quadrature event does not open new physical IR phase space. Angular endpoint resolution is an independent, dominant numerical uncertainty. Increase nmu or use a verified endpoint-adapted quadrature and recheck all moments.",
    }


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--quick",action="store_true")
    parser.add_argument("--output",default="outputs/angular_endpoint_audit.json")
    opts=parser.parse_args()
    result=run(quick=opts.quick)
    p=Path(opts.output);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(result,indent=2),encoding="utf-8")
    print(json.dumps(result,indent=2))
