"""Independent initial-vacuum source quadrature audit of adaptive closed kinetics.

The continuum on-shell EMPTY-daughter parent source can be integrated
exactly in 1D: g1*M/(2*pi^2)*int_0^inf log(1+exp[-q-M^2/(4q)]) dq.
For finite-cell closed triplets with shared mother energy E_i+E_j,
the initial decay rate is sum(K_ij * fFD(E_i+E_j)).
Compare it BOTH to the continuum 1D integral and to the actual ODE
d[n_phi]/d(t Gamma) at t=0.

The finite-cell rate can converge non-monotonically as the daughter
bins/centroid-on-shell mask change. Correct source normalization is not
sufficient for nonlinear/stimulated f(q) convergence.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import expit

from adaptive_resonant_lattice import AdaptiveClosedDecay
from resonant_graph_repair import bridged_indices
from onshell_collision_2d import vacuum_reference
from schur_closed_bdf import solve_closed_schur
from sparse_kinetic_jacobian import SparseClosedDecay


def one_source(core=100,stride=50,max_index=8500,mass=.1):
    indices=bridged_indices(core=core,stride=stride,max_index=max_index)
    model=AdaptiveClosedDecay(indices,M1_T=mass)
    n=model.nq
    fparent=expit(-model.parent_E)
    y=np.r_[np.zeros(2*n),fparent]
    rate_by_reaction=float(np.sum(model.K*fparent[model.parent_index]))
    dydt=model.rhs(0.,y)
    rate_by_ode=float(np.dot(model.vol_phi,dydt[:n]))
    exact=vacuum_reference(M1_T=mass)
    if abs(rate_by_reaction-rate_by_ode)>1e-11:
        raise ArithmeticError("Closed ODE initial source disagrees with its collision integral")
    return {
      "mesh":{"core":core,"stride":stride,"max_index":max_index,
              "delta":model.delta,"M1_over_T":mass},
      "daughter_modes":n,"parent_modes":model.nparent,
      "daughter_domain_upper_edge":float(model.edges[-1]),
      "reaction_graph_invariants":model.graph["active_graph_collision_invariants"],
      "initial_grid_source_rate":rate_by_reaction,
      "ODE_initial_phi_number_rate":rate_by_ode,
      "continuum_analytic_source_rate":exact,
      "grid_source_relative_error":float(rate_by_reaction/exact-1),
      "initial_parent_number_relative_error":
         float(model.mother_inventory_relative_bias),
      "warning":"Vacuum/empty daughters only; not the non-linear stimulated integrated production."
    }


def run(quick=False,dynamic=False):
    qmax=(5000,7000) if quick else (5000,7000,8500,10000)
    rows=[one_source(max_index=m) for m in qmax]
    alternative_meshes=[] if quick else [
         one_source(core=100,stride=25,max_index=10000),
         one_source(core=120,stride=60,max_index=10000),
         one_source(core=150,stride=50,max_index=10000)]
    dynamic_rows=[]
    if dynamic:
        for m in qmax:
            model=SparseClosedDecay(bridged_indices(max_index=m))
            out=solve_closed_schur(model,tend=.001)
            dynamic_rows.append({
                "max_index":m,"short_time_tGamma":.001,
                "total_phi_over_T3":out["n_phi_over_T3"],
                "mean_phi_q":out["mean_phi_q"],
                "not_frozen":True
            })
    sensitivity=None
    if dynamic and len(rows)>=2:
        initial_change=rows[-1]["initial_grid_source_rate"]/rows[-2]["initial_grid_source_rate"]-1
        evolved_change=dynamic_rows[-1]["total_phi_over_T3"]/dynamic_rows[-2]["total_phi_over_T3"]-1
        sensitivity={
          "last_two_max_indices":[qmax[-2],qmax[-1]],
          "initial_vacuum_rate_fractional_change":float(initial_change),
          "nonlinear_tGamma0p001_phi_number_fractional_change":float(evolved_change),
          "ratio_of_two_cutoff_sensitivities":float(
              abs(evolved_change)/max(abs(initial_change),1e-15)),
          "interpretation":"Numerical sensitivity ratio only. Not a physical Bose enhancement factor; the nonlinear histories, source depletion and inverse reactions differ."
        }
    return {
       "title":"Continuum-vs-discrete on-shell EMPTY-daughter vacuum source",
       "physical_model":"M1/T=.1 massless daughters thermal FD initial chi1, g1=4",
       "source_quad_method":"analytic exact hyperbolic cell areas + centroid on-shell mask + unique mother state; compared to 1D vacuum reference",
       "results":rows,"alternative_mesh_topologies":alternative_meshes,
       "optional_short_time_closed_dynamics":dynamic_rows,
       "optional_cutoff_sensitivity_comparison":sensitivity,
       "caveat":"A small initial source quadrature error does not bound the high-occupation nonlinear collision dynamics or imply a physical FRW relic PSD.",
    }


if __name__=="__main__":
    pa=argparse.ArgumentParser()
    pa.add_argument("--quick",action="store_true")
    pa.add_argument("--dynamic",action="store_true")
    pa.add_argument("--output",default="outputs/vacuum_source_convergence.json")
    arg=pa.parse_args()
    ans=run(arg.quick,arg.dynamic)
    dst=Path(arg.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(ans,indent=2),encoding="utf-8")
    print(json.dumps(ans,indent=2))
