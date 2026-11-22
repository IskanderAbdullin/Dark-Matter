"""Self-convergence dt scan for nonthermal quantum Zprime nu+antinu Strang.

For a FIXED finite collision graph, compare max absolute mode differences
||Y_dt - Y_dt/2||_infty for all three species. Ratios ~4 across halving
mean second temporal order; they DO NOT imply continuum quadrature or a
cosmological physical rate. No thermal reset, Γ0=1 synthetic.
"""
from __future__ import annotations
import argparse
import json
import math
import numpy as np
from kinetics.vector_neutrino_exact_event_strang import pilot


def scan(h=.1,tau=1.,dts=(.125,.0625,.03125,.015625,.0078125,
                           .00390625,.001953125,.0009765625)):
    if not dts or any(x<=0 for x in dts) or any(
            abs(round(tau/x)*x-tau)>1e-10 for x in dts):
        raise ValueError("positive commensurate dt")
    states={}
    rows=[]
    for dt in dts:
        (b,fn,fa),out=pilot(h=h,total_tau=tau,dt=dt)
        states[dt]=(b,fn,fa)
        rows.append({
           "dt":dt,"steps":out["steps"],
           "low_vector_b":out["low_vector_b"],
           "min_step_entropy_gain":out["min_step_entropy_gain"],
           "max_abs_charge_or_energy_residual":max(out["max_C_drift"],
                   out["max_L_drift"],out["max_E_drift"]),
           "min_b":out["min_b"],"min_fermion":out["min_fermion"],
           "max_fermion":out["max_fermion"]})
    differences=[]
    prev=None
    for dt in dts[:-1]:
        fine=dt/2
        if fine not in states:raise ValueError("dts must be adjacent halving")
        diff=max(float(np.max(np.abs(a-b))) for a,b in
                 zip(states[dt],states[fine]))
        differences.append({"dt":dt,"error_to_half_dt":diff,
          "successive_error_ratio":None if prev is None else prev/diff})
        prev=diff
    return {
      "finite_graph_only":True,"h":h,"tau_in_Gamma0_inv":tau,
      "mV":.4,"Emax":4,
      "rows":rows,"step_halving_differences":differences,
      "not_physical_chi_gamma_el_or_86eV_relic":True,
      "no_cosmological_FRW_or_continuum_claim":True}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--h",type=float,default=.1)
    parser.add_argument("--tau",type=float,default=1.)
    args=parser.parse_args()
    print(json.dumps(scan(h=args.h,tau=args.tau),indent=2))
