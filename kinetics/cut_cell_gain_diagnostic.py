"""Decompose on-shell finite-cell vacuum-source error by exact cut-cell FD quadrature.

The continuum source in a daughter bin [xlo,xhi] x [ylo,yhi],
restricted to xy>=c=M1²/4, is

 C int_dx [ log(1+exp[-x-max(ylo,c/x)])
            -log(1+exp[-x-yhi]) ] , C=g1*M1/(2*pi²).

Split x into a partially cut interval (c/yhi,c/ylo) and a fully
kinematically allowed interval (c/ylo,xhi). Fixed Gauss-Legendre nodes
evaluate a smooth 1D integrand on each segment.

This is a DIAGNOSTIC only: it does not change any collision channel,
mother occupancy, physical parent energy, H theorem, or dynamical ODE.
In particular, cells with a forbidden centroid cannot simply be
included as resonant parent energies without breaking the on-shell
energy relation. They are reported as missing source separately.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import expit,roots_legendre

from adaptive_resonant_lattice import AdaptiveClosedDecay
from resonant_graph_repair import bridged_indices
from uv_source_decomposition import exact_rectangle_source


def integrate_cut_cells(edges, mass=.1, quadrature_n=16):
    edges=np.asarray(edges,dtype=float)
    if len(edges)<3 or edges[0]<0 or np.any(np.diff(edges)<=0):
        raise ValueError("Need sorted nonnegative daughter-energy cell edges")
    if mass<=0 or quadrature_n<2:
        raise ValueError("Positive mass and Gauss nodes >=2 required")
    lo=edges[:-1]
    hi=edges[1:]
    N=len(lo)
    xl=np.broadcast_to(lo[:,None],(N,N)).ravel()
    xh=np.broadcast_to(hi[:,None],(N,N)).ravel()
    yl=np.broadcast_to(lo[None,:],(N,N)).ravel()
    yh=np.broadcast_to(hi[None,:],(N,N)).ravel()
    c=mass*mass/4
    cutoff_ylo=np.divide(c,yl,out=np.full_like(yl,np.inf),where=yl>0)
    first=np.maximum(xl,c/yh)
    last=np.minimum(xh,cutoff_ylo)
    start_full=np.maximum(xl,cutoff_ylo)
    # Each interval is zero unless right edge is strictly larger;
    # reset invalid endpoints to safe finite numbers before evaluation.
    partial_valid=last>first
    full_valid=xh>start_full
    pa=np.where(partial_valid,first,xl)
    pb=np.where(partial_valid,last,xl)
    fa=np.where(full_valid,start_full,xl)
    fb=np.where(full_valid,xh,xl)
    nodes,weights=roots_legendre(quadrature_n)

    def integrate_piece(left,right,partial):
        length=right-left
        X=(left[:,None]+.5*length[:,None]*(1+nodes[None,:]))
        Ylo=c/np.maximum(X,1e-300) if partial else yl[:,None]
        inner=np.logaddexp(0.,-X-Ylo)-np.logaddexp(0.,-X-yh[:,None])
        result=.5*length*np.sum(inner*weights[None,:],axis=1)
        return result

    integrated=integrate_piece(pa,pb,True)+integrate_piece(fa,fb,False)
    integrated=np.maximum(integrated,0.)
    return integrated.reshape(N,N)


def evaluate_cells(max_index=8500,core=100,stride=50,
                   mass=.1,quadrature_n=16):
    model=AdaptiveClosedDecay(
        bridged_indices(core=core,stride=stride,max_index=max_index),
        M1_T=mass)
    source_factor=4*mass/(2*math.pi**2)
    exact_cell=source_factor*integrate_cut_cells(
        model.edges,mass=mass,quadrature_n=quadrature_n)
    active=model.K>0
    true_on_domain=exact_rectangle_source(float(model.edges[-1]),mass=mass)
    total_cell=float(exact_cell.sum())
    included=float(np.sum(exact_cell[active]))
    excluded=float(np.sum(exact_cell[~active]))
    sampled=float(np.sum(model.K*expit(
        -(model.e[:,None]+model.e[None,:]))))
    area_integral_relative_error=total_cell/true_on_domain-1
    quadrature_error_over_domain=(sampled-included)/true_on_domain
    centroid_mask_missing_fraction=excluded/true_on_domain
    total_error=sampled/true_on_domain-1
    closure=quadrature_error_over_domain-centroid_mask_missing_fraction-total_error
    assert abs(closure)<1e-8,closure
    if abs(area_integral_relative_error)>2e-6:
        raise AssertionError("Direct cut-cell sum does not reproduce rectangle continuum integral")
    return {
       "mesh":{"core":core,"stride":stride,"max_index":max_index,"mass":mass},
       "daughter_modes":model.nq,
       "active_reaction_pairs":model.active_pairs,
       "quadrature_n":quadrature_n,
       "topology_active_collision_invariants":
           model.graph["active_graph_collision_invariants"],
       "exact_rectangle_source":true_on_domain,
       "all_cut_cell_integral":total_cell,
       "relative_cut_cell_quadrature_error":area_integral_relative_error,
       "kinematically_valid_centroid_cell_integral":included,
       "kinematically_forbidden_centroid_cell_integral":excluded,
       "sampled_kinetic_vacuum_source":sampled,
       "relative_error_due_to_centroid_fFD_approximation":
           float(quadrature_error_over_domain),
       "relative_error_due_to_omitted_cut_cells":
           float(-centroid_mask_missing_fraction),
       "relative_total_finite_cell_source_error":float(total_error),
       "error_decomposition_residual":float(closure),
       "warning":"Continuum thermal VACUUM source diagnostic only. Correcting the masked reaction network requires preserving parent on-shell energy, three invariants and quantum entropy, not just adding missing weights."
    }


def run(quick=False):
    qmax=(5000,) if quick else (5000,7000,8500,10000)
    rows=[evaluate_cells(max_index=x) for x in qmax]
    accuracy=(
        evaluate_cells(max_index=5000,quadrature_n=32)
        if quick else evaluate_cells(max_index=8500,quadrature_n=32))
    comparison=next(v for v in rows if v["mesh"]["max_index"]==
                    accuracy["mesh"]["max_index"])
    relative_sensitivity=abs(accuracy["all_cut_cell_integral"]/
                             comparison["all_cut_cell_integral"]-1)
    assert relative_sensitivity<1e-6
    return {
       "scope":"Source-only FD integration of all on-shell cut cells, including discarded forbidden-centroid channels",
       "physical_params":"M1/T=.1, massless daughters, initial thermal FD mother, empty daughter states",
       "cases":rows,
       "higher_quadrature_test":accuracy,
       "relative_change_gauss16_to_gauss32":relative_sensitivity,
       "no_collision_operator_change":True,
       "remaining":"Design a physically resonant and entropic cut-cell discretization; full nonlinear and FRW convergence not established."
    }


if __name__=="__main__":
    pa=argparse.ArgumentParser()
    pa.add_argument("--quick",action="store_true")
    pa.add_argument("--output",default="outputs/cut_cell_vacuum_diagnostic.json")
    arg=pa.parse_args()
    result=run(quick=arg.quick)
    destination=Path(arg.output)
    destination.parent.mkdir(parents=True,exist_ok=True)
    destination.write_text(json.dumps(result,indent=2),encoding="utf-8")
    print(json.dumps(result,indent=2))
