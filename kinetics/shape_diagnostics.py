"""Number-weighted frozen-spectrum SHAPE diagnostics for finite-bin Bose kinetics.

These diagnostics analyze a finite-bin snapshot; they do not make it
a frozen cosmological relic. CDF/quantiles integrate the shell's q^2
phase-space measure under the piecewise-constant-bin approximation.

A high Bose occupation in a tiny q bin is grid-dependent; do not infer
a zero-mode condensate or physical CLASS initial distribution from it.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np

from sparse_kinetic_jacobian import SparseClosedDecay
from schur_closed_bdf import solve_closed_schur
from resonant_graph_repair import bridged_indices


CHECK_Q=(0.0005,0.001,0.002,0.003,0.005,0.01,0.03,0.1)


def validate_snapshot(spectrum):
    edges=np.asarray(spectrum["q_bin_edges"],dtype=float)
    b=np.asarray(spectrum["phi_occupations"],dtype=float)
    n=np.asarray(spectrum["phi_number_per_bin"],dtype=float)
    if len(edges)!=len(b)+1 or len(n)!=len(b):
        raise ValueError("Inconsistent snapshot bin dimensions")
    if np.any(edges<0) or np.any(np.diff(edges)<=0):
        raise ValueError("q-bin edges must be strictly increasing and nonnegative")
    if np.any(b<0) or np.any(n<0) or not np.all(np.isfinite(n)):
        raise ValueError("occupations/number bins must be nonnegative finite")
    total=float(sum(n))
    if total<=0:
        raise ValueError("Cannot normalize an empty Bose population")
    if abs(total/spectrum["total_phi_number"]-1)>1e-8:
        raise ArithmeticError("Bin number moment inconsistent with kinetic source")
    return edges,b,n,total


def cdf_number_at_q(edges,n_bin,q):
    edges=np.asarray(edges,dtype=float)
    n_bin=np.asarray(n_bin,dtype=float)
    if q<0:
        raise ValueError("momentum must be nonnegative")
    left=edges[:-1];right=edges[1:]
    shell_fraction=np.clip((q**3-left**3)/(right**3-left**3),0,1)
    return float(np.dot(n_bin,shell_fraction)/np.sum(n_bin))


def quantile_number_q(edges,n_bin,probability):
    edges=np.asarray(edges,dtype=float)
    n_bin=np.asarray(n_bin,dtype=float)
    if not 0<=probability<=1:
        raise ValueError("quantile probability outside [0,1]")
    total=float(n_bin.sum())
    if total<=0:
        raise ValueError("empty number spectrum")
    if probability==0:
        return float(edges[np.flatnonzero(n_bin>0)[0]])
    if probability==1:
        return float(edges[np.flatnonzero(n_bin>0)[-1]+1])
    cumul=np.cumsum(n_bin)
    index=min(int(np.searchsorted(cumul,probability*total)),len(n_bin)-1)
    while n_bin[index]<=0 and index+1<len(n_bin):
        index+=1
    before=cumul[index-1] if index else 0.
    fraction=np.clip((probability*total-before)/n_bin[index],0,1)
    return float((edges[index]**3+fraction*(
        edges[index+1]**3-edges[index]**3))**(1/3))


def describe_number_spectrum(spectrum,check_q=CHECK_Q):
    edges,b,number,total=validate_snapshot(spectrum)
    thresholds={f"{q:.7g}":cdf_number_at_q(edges,number,q) for q in check_q}
    quantiles={f"q{round(p*100):02d}":quantile_number_q(edges,number,p)
               for p in (.01,.1,.5,.9,.99)}
    return {
        "number_over_T3":total,
        "quantile_momenta":quantiles,
        "normalized_number_CDF":thresholds,
        "fraction_q_below_0p1":cdf_number_at_q(edges,number,.1),
        "highest_bin_occupation":float(max(b)),
        "bin_count":len(b),
        "warning":"Piecewise-constant shell distribution at short tGamma, not primordial freezeout or evidence of condensation."
    }


def one_case(max_index,tend=.001):
    model=SparseClosedDecay(bridged_indices(max_index=max_index))
    answer=solve_closed_schur(model,tend=tend,return_spectrum=True)
    report=describe_number_spectrum(answer["daughter_spectrum"])
    assert abs(report["number_over_T3"]/answer["n_phi_over_T3"]-1)<1e-12
    assert abs(report["fraction_q_below_0p1"]-
               answer["phi_fraction_below_q0p1"])<1e-8
    return {"upper_daughter_index":max_index,"summary":report,
            "number_and_energy_residuals":answer["number_and_energy_residuals"],
            "short_time_only":True}


def run(quick=False):
    cutoffs=(5000,) if quick else (5000,7000,8500,10000)
    rows=[one_case(m) for m in cutoffs]
    shape_delta=None
    if len(rows)>=2:
        a,b=rows[-2]["summary"],rows[-1]["summary"]
        keys=a["normalized_number_CDF"].keys()
        shape_delta={
            "sup_CDF_difference_over_reported_fixed_thresholds":
                max(abs(a["normalized_number_CDF"][k]-
                        b["normalized_number_CDF"][k]) for k in keys),
            "q50_relative_difference":float(
                b["quantile_momenta"]["q50"]/a["quantile_momenta"]["q50"]-1),
            "q90_relative_difference":float(
                b["quantile_momenta"]["q90"]/a["quantile_momenta"]["q90"]-1),
            "cutoffs_compared":[rows[-2]["upper_daughter_index"],
                                rows[-1]["upper_daughter_index"]],
        }
    return {
        "status":"FINITE-BIN SHORT-TIME NORMALIZED SHAPE CHECK — NOT COSMOLOGICAL PSD",
        "time_tGamma":.001,"statistics":"q² dq number measure within each finite shell",
        "reference":"closed 3-invariant 1<->2 mother-depleted Schur BDF",
        "rows":rows,"shape_difference":shape_delta,
        "prohibition":"No interpretation as primordial 86-eV DM spectrum, no BEC claim, do not pass to CLASS.",
    }


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--quick",action="store_true")
    p.add_argument("--output",default="outputs/closed_short_time_shape.json")
    opts=p.parse_args()
    out=run(opts.quick)
    path=Path(opts.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,indent=2),encoding="utf-8")
    print(json.dumps(out,indent=2))
