"""Separate finite daughter-energy domain truncation from quadrature error.

For massless daughters, M1/T=m, c=m²/4, the exact 1<->2
EMPTY-daughter on-shell source in rectangular 0<q_phi,q_chi<Q is:

 R_rect(Q) = g1*m/(2*pi²) * integral_(c/Q)^Q dq [
     log(1+exp(-q-c/q)) - log(1+exp(-q-Q)) ].

This removes the common mistake of attributing all UV dependence to
nonlinear integration or finite-cell accuracy. It is an independent
continuum truncated-domain calculation; no adaptive binning used.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad

from onshell_collision_2d import vacuum_reference
from vacuum_source_convergence import one_source


def exact_rectangle_source(Q, mass=.1, g1=4):
    """Exactly integrate χ2 analytically, φ with 1D adaptive quadrature."""
    if Q<=0 or mass<=0 or g1<=0:
        raise ValueError("Q,mass,g1 must be positive")
    c=mass*mass/4
    if Q*Q<=c:
        return 0.0
    lower=c/Q
    def integrand(q):
        a=float(np.logaddexp(0.,-q-c/q))
        b=float(np.logaddexp(0.,-q-Q))
        return a-b
    integral=quad(integrand,lower,Q,epsabs=1e-13,
                  epsrel=1e-11,limit=300)[0]
    return g1*mass/(2*math.pi**2)*integral


def one_comparison(max_index,core=100,stride=50,mass=.1):
    source=one_source(core=core,stride=stride,
                      max_index=max_index,mass=mass)
    Q=source["daughter_domain_upper_edge"]
    full=source["continuum_analytic_source_rate"]
    truncated=exact_rectangle_source(Q,mass=mass)
    discrete=source["initial_grid_source_rate"]
    cutoff=truncated/full-1.
    cell=discrete/truncated-1.
    total=discrete/full-1.
    multiplicative_residual=(1+cutoff)*(1+cell)-(1+total)
    if abs(multiplicative_residual)>1e-12:
        raise ArithmeticError("UV vs finite-cell error decomposition failed")
    return {
        "mesh":source["mesh"],
        "daughter_energy_max_EDGE":Q,
        "daughter_modes":source["daughter_modes"],
        "continuum_infinite_domain_source":full,
        "continuum_finite_rectangle_source":truncated,
        "finite_cell_on_shell_source":discrete,
        "relative_domain_truncation_error":float(cutoff),
        "relative_cell_rate_error_within_same_domain":float(cell),
        "relative_total_source_error":float(total),
        "exact_multiplicative_decomposition_residual":float(multiplicative_residual),
        "initial_mother_FD_number_relative_error":
            source["initial_parent_number_relative_error"],
        "not_a_nonlinear_source_convergence_test":True
    }


def run(quick=False):
    indices=(5000,8500) if quick else (5000,7000,8500,10000)
    rows=[one_comparison(q) for q in indices]
    reference=vacuum_reference()
    assert all(abs(r["continuum_infinite_domain_source"]/reference-1)<1e-12
               for r in rows)
    assert all(r["relative_domain_truncation_error"]<0 for r in rows)
    assert all(abs(r["exact_multiplicative_decomposition_residual"])<1e-12
               for r in rows)
    alternative=[] if quick else [
        one_comparison(10000,core=100,stride=25),
        one_comparison(10000,core=120,stride=60),
        one_comparison(10000,core=150,stride=50)]
    return {
      "model":"on-shell 1<->2 massless daughters, FD mother, initially empty daughters",
      "physical_scope":"exact CONTINUUM rectangular-domain vacuum rate vs finite-cell initial rate",
      "constant_temperature_M1_over_T":.1,
      "full_source_continuum_reference":reference,
      "cases":rows,
      "alternative_meshing_at_same_approximate_Q":alternative,
      "main_result":"Dominant low-Q error is excluded daughter momentum phase space. At high Q a residual finite-cell error remains; neither is a bound for nonlinear stimulated dynamics.",
      "caveat":"Not a frozen primordial PSD, dynamic Bose gain or physically self-consistent FRW scenario."
    }


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--quick",action="store_true")
    p.add_argument("--output",default="outputs/uv_source_error_decomposition.json")
    opts=p.parse_args()
    result=run(opts.quick)
    destination=Path(opts.output)
    destination.parent.mkdir(parents=True,exist_ok=True)
    destination.write_text(json.dumps(result,indent=2),encoding="utf-8")
    print(json.dumps(result,indent=2))
