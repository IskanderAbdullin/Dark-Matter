"""Positive thermally calibrated rate on the SAME on-shell resonant 1<->2 graph.

For active centroid-on-shell reaction i+j<->k with E_k=e_i+e_j,
replace K_ij= C * area_ij by

  Kcal_ij = C * [integral over physically allowed PART of cell
                  fFD(x+y) dxdy] / fFD(E_k).

This exactly integrates the INITIAL thermal empty-daughter forward
source ON THE RETAINED REACTION GRAPH. Cells whose CENTROID is
kinematically forbidden remain excluded: we do not manufacture
unphysical resonances or duplicate parent states to fix totals.

Kcal>=0, and BOTH forward/reverse occupation factors use the SAME
Kcal. Therefore every discrete reaction retains exact energy,
two number invariants, quantum boundary positivity, thermal detailed
balance and H theorem. This is not an exact continuum collision kernel
for arbitrary nonequilibrium distributions (reference-FD calibration).
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import expit

from cut_cell_gain_diagnostic import integrate_cut_cells
from closed_three_species_lattice import allowed_area_from_zero
from sparse_kinetic_jacobian import SparseClosedDecay
from resonant_graph_repair import bridged_indices
from schur_closed_bdf import solve_closed_schur


class ThermalCutCellCalibratedDecay(SparseClosedDecay):
    def __init__(self,indices,delta=.001,M1_T=.1,g1=4,gchi=4,
                 quadrature_n=16):
        super().__init__(indices,delta=delta,M1_T=M1_T,
                         g1=g1,gchi=gchi)
        old_mask=self.K>0
        full_thermal_source=integrate_cut_cells(
            self.edges,mass=M1_T,quadrature_n=quadrature_n)
        bath=expit(-(self.e[:,None]+self.e[None,:]))
        if np.any(bath[old_mask]<=0):
            raise ArithmeticError("Reference bath underflow")
        thermal_geometric=full_thermal_source[old_mask]
        if np.any(thermal_geometric<=0):
            raise ArithmeticError("Positive on-shell bin missing thermal source")
        C=g1*M1_T/(2*np.pi**2)
        new_K=np.zeros_like(self.K)
        new_K[old_mask]=C*thermal_geometric/bath[old_mask]
        if np.any(new_K[old_mask]<=0):
            raise ArithmeticError("Calibration must preserve positive collisions")
        if np.any((new_K>0)!=old_mask):
            raise ArithmeticError("Calibration changed collision graph")
        self.K=new_K
        self.w=self.K[self.i,self.j]
        self.reference_quadrature_n=quadrature_n
        self.initial_source_calibrated=float(np.sum(self.K*bath))
        self.reference_thermal_source_on_allowed_centroids=float(
            C*full_thermal_source[old_mask].sum())
        self.excluded_centroid_cutcell_source=float(
            C*full_thermal_source[~old_mask].sum())
        self.unmodified_source=float(np.sum(
            (g1*M1_T/(2*np.pi**2)) *
            allowed_area_from_zero(self.edges,M1_T)*old_mask*bath))
        self.surrogate_not_class_ready=True

    def audit(self):
        relative=abs(self.initial_source_calibrated/
                     self.reference_thermal_source_on_allowed_centroids-1)
        if relative>1e-12:
            raise AssertionError("Thermal cut-cell calibration failed")
        structural=self.thermodynamic_checks()
        return {
           "reference_quadrature_nodes":self.reference_quadrature_n,
           "daughter_count":self.nq,
           "shared_parent_modes":self.nparent,
           "retained_reaction_count":self.active_pairs,
           "graph_invariants":self.graph["active_graph_collision_invariants"],
           "old_thermally_empty_source":self.unmodified_source,
           "new_calibrated_empty_source":self.initial_source_calibrated,
           "continuum_source_on_center_allowed_cells":
               self.reference_thermal_source_on_allowed_centroids,
           "continuum_source_omitted_because_centroid_offshell":
               self.excluded_centroid_cutcell_source,
           "calibration_to_allowed_cells_relative_error":relative,
           "thermal_detailed_balance_and_h_checks":structural,
           "warning":"No physical guarantee for nonlinear off-thermal Bose spectrum. Off-shell centroid cut cells are still excluded, and full FRW/zero-mode physics absent."
        }


def run(quick=False):
    cutoffs=(5000,) if quick else (5000,8500,10000)
    results=[]
    for end in cutoffs:
        indices=bridged_indices(max_index=end)
        m=ThermalCutCellCalibratedDecay(indices)
        audit=m.audit()
        conventional=SparseClosedDecay(indices)
        if quick:
            time=.0002
        else:
            time=.001
        old=solve_closed_schur(conventional,tend=time)
        calibrated=solve_closed_schur(m,tend=time)
        results.append({
            "max_index":end,
            "audit":audit,
            "t_Gamma":time,
            "uncalibrated_nphi_over_T3":old["n_phi_over_T3"],
            "thermally_calibrated_nphi_over_T3":
                calibrated["n_phi_over_T3"],
            "relative_dynamic_nphi_change":float(
                calibrated["n_phi_over_T3"]/old["n_phi_over_T3"]-1),
            "uncalibrated_mean_phi_q":old["mean_phi_q"],
            "calibrated_mean_phi_q":calibrated["mean_phi_q"],
            "calibrated_charge_energy_residuals":
                calibrated["number_and_energy_residuals"]
        })
    return {
        "status":"POSITIVE CUT-CELL THERMAL CALIBRATION; source-only normalization correction in a CLOSED finite-grid dynamical test",
        "physics_scope":"M1/T=.1, reference FD mother at T=1, empty daughters initially, Minkowski; exact resonant representative energies",
        "graph_changed":False,
        "physical_equations_added":False,
        "not_a_frozen_dark_matter_PSD":True,
        "cases":results,
        "caveat":"Calibration is exactly correct only for empty-daughter thermal source on active centroid-resonant cells. It cannot capture omitted on-shell phase space without a new energy-resonant discretization or prove nonlinear IR/FRW convergence.",
    }


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--quick",action="store_true")
    p.add_argument("--output",default="outputs/thermal_cutcell_calibration.json")
    opt=p.parse_args()
    d=run(quick=opt.quick)
    path=Path(opt.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(d,indent=2),encoding="utf-8")
    print(json.dumps(d,indent=2))
