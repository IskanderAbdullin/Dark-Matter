"""Small-X, topology-rebuilt collision-plus-Liouville FRW SPLIT PILOT.

Advance a CLOSED mother/daughter 1<->2 system between x=M/T values
on a prescribed radiation background with constant entropy degrees.
At every x, rebuild the on-shell reaction graph K(x), ensuring each
collisional mother energy E_k=E_phi_i+E_chi_j exactly. A mother state
is reinterpreted on a MOVING comoving momentum grid and conservatively
rebinned in q**3 between x values. Daughters (massless) have fixed q=p/T.

This is a deliberately NARROW structural experiment. The shared
parent resonant-energy centers must remain above x, and the set of
parent energies must stay unchanged, or fail closed. The collision
step conserves two charges and discrete energy; remapping conserves
both charges/Pauli but DOES NOT conserve physical energy exactly.
Physical free-redshift work, collision discretization error and
remapping energy error are tracked SEPARATELY. Not a cosmological
relic-density, zero-mode, or CLASS calculation.

Radiation H(x)=H(x=1)/x², d(tGamma)/dx=kappa*x.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expit

from sparse_kinetic_jacobian import SparseClosedDecay
from resonant_graph_repair import bridged_indices
from parent_q_remap import (
    moving_parent_edges, thermal_shell_averages, rebin, number, energy
)


class NarrowClosedFRWPilot:
    def __init__(self,x0=.02,xend=.0202,steps=4,kappa=500.,
                 daughter_max_index=300,allow_mother_mode_changes=False):
        if x0<=0 or xend<=x0 or steps<1 or kappa<=0:
            raise ValueError("Positive narrow FRW range, steps, coupling required")
        self.xgrid=np.linspace(x0,xend,steps+1)
        self.kappa=kappa
        self.allow_mother_mode_changes=bool(allow_mother_mode_changes)
        self.indices=bridged_indices(max_index=daughter_max_index)
        try:
            self.models=[SparseClosedDecay(self.indices,M1_T=float(x))
                         for x in self.xgrid]
        except AssertionError as exc:
            if "Spurious invariant graph" not in str(exc):
                raise
            raise ValueError("Extra or missing reaction invariants: refuse evolution") from exc
        reference=self.models[0]
        for model,xx in zip(self.models,self.xgrid):
            if not self.allow_mother_mode_changes and not np.array_equal(model.parent_E,reference.parent_E):
                raise ValueError("Parent resonant energies changed: explicit graph rebuild required")
            if model.graph["active_graph_collision_invariants"]!=3:
                raise ValueError("Extra or missing reaction invariants: refuse evolution")
            if xx>=model.parent_E[0]:
                raise ValueError("Parent energy center below instantaneous physical rest mass")
        # One global COMOVING momentum domain for every moving grid.
        # Parent energy modes may appear/disappear if and only if their
        # new exact on-shell graph independently passes the rank gate.
        self.qmax=max(math.sqrt((model.parent_E[-1]+.5*(model.parent_E[-1]-model.parent_E[-2]))**2-xx**2)
                      for model,xx in zip(self.models,self.xgrid))+1e-12
        self.parent_edges=[]
        for model,x in zip(self.models,self.xgrid):
            qedge=moving_parent_edges(model.parent_E,float(x),self.qmax)
            self.parent_edges.append(qedge)
            # Make parent number volumes match the SAME global comoving
            # momentum domain at every x, rather than inheriting an
            # x-dependent last-bin cutoff of the flat-space model.
            model.vol_parent=(4/(6*np.pi**2))*np.diff(qedge**3)
        self.n=reference.nq

    @staticmethod
    def _charges(model,state):
        m=model.moments(state)
        return (m["charge_chi1_plus_phi"],m["charge_chi1_plus_chi2"])

    def physical_rho(self,model,parent_edges,state,x):
        b,f,g=model.unpack(state)
        daughter=float(np.dot(model.vol_phi*model.e,b)+
                       np.dot(model.vol_chi*model.e,f))
        mother=energy(parent_edges,g,x)
        return daughter+mother

    def _collide(self,model,state,dtau,rtol=2e-7,atol=1e-12):
        result=solve_ivp(model.rhs_sparse,(0,dtau),state,method="BDF",
                         jac=model.jacobian,rtol=rtol,atol=atol,
                         max_step=max(dtau/8,1e-12))
        if not result.success:
            raise RuntimeError(result.message)
        last=result.y[:,-1]
        b,f,g=model.unpack(last)
        if np.min(last)<-1e-8 or np.max(f)>1+1e-8 or np.max(g)>1+1e-8:
            raise ArithmeticError("Quantum occupation bounds broken")
        return last,{"nfev":int(result.nfev),"njev":int(result.njev),"nlu":int(result.nlu)}

    def run(self):
        initial_model=self.models[0]
        initial_g=thermal_shell_averages(self.parent_edges[0],float(self.xgrid[0]))
        state=np.r_[np.zeros(2*self.n),initial_g]
        first_charges=self._charges(initial_model,state)
        rows=[]
        for step in range(len(self.models)-1):
            m=self.models[step]
            nxt=self.models[step+1]
            x=float(self.xgrid[step]);xnew=float(self.xgrid[step+1])
            before=self._charges(m,state)
            before_discrete=m.moments(state)["total_discrete_energy"]
            phys_before=self.physical_rho(m,self.parent_edges[step],state,x)
            dtau=self.kappa*(xnew*xnew-x*x)/2
            after_collision,counts=self._collide(m,state,dtau)
            collision_discrete_error=(
                m.moments(after_collision)["total_discrete_energy"]-before_discrete)
            collision_physical_error=(
                self.physical_rho(m,self.parent_edges[step],after_collision,x)-
                phys_before)
            b,f,g=m.unpack(after_collision)
            phys_free_redshift=(
                self.physical_rho(m,self.parent_edges[step],after_collision,xnew)-
                self.physical_rho(m,self.parent_edges[step],after_collision,x))
            remapped=rebin(self.parent_edges[step],
                           self.parent_edges[step+1],g)
            state=np.r_[b.copy(),f.copy(),remapped]
            after=self._charges(nxt,state)
            remap_energy_error=(
                self.physical_rho(nxt,self.parent_edges[step+1],state,xnew)-
                self.physical_rho(m,self.parent_edges[step],after_collision,xnew))
            charge_drift=max(abs(after[0]-before[0]),abs(after[1]-before[1]))
            if abs(collision_discrete_error)>1e-8 or charge_drift>1e-8:
                raise ArithmeticError("Collision/remap charge or discrete energy violated")
            rows.append({
                "x_start":x,"x_end":xnew,"delta_tGamma":dtau,
                "active_on_shell_channels_start":m.active_pairs,
                "active_on_shell_channels_end":nxt.active_pairs,
                "parent_states":m.nparent,
                "next_parent_states":nxt.nparent,
                "parent_mode_count_change":nxt.nparent-m.nparent,
                "before_charge1":before[0],"before_charge2":before[1],
                "after_charge1":after[0],"after_charge2":after[1],
                "maximum_charge_drift":charge_drift,
                "collision_discrete_energy_error":collision_discrete_error,
                "collision_physical_shell_energy_error":collision_physical_error,
                "collisionless_FRW_energy_change":phys_free_redshift,
                "finite_volume_remap_energy_error":remap_energy_error,
                "min_fermion_occupation":float(min(f.min(),remapped.min())),
                "max_fermion_occupation":float(max(f.max(),remapped.max())),
                **counts
            })
        final=self.models[-1].moments(state)
        max_initial_drift=max(abs(final["charge_chi1_plus_phi"]-first_charges[0]),
                              abs(final["charge_chi1_plus_chi2"]-first_charges[1]))
        if max_initial_drift>1e-8:raise ArithmeticError("Global closed charges drifted")
        return {
            "status":"NARROW SPLIT-STEP CONSTANT-GSTAR FRW STRUCTURAL BENCHMARK",
            "x0":float(self.xgrid[0]),"xend":float(self.xgrid[-1]),
            "kappa":self.kappa,"daughter_modes":self.n,
            "shared_parent_modes":self.models[0].nparent,
            "final_shared_parent_modes":self.models[-1].nparent,
            "parent_mode_changes_permitted":self.allow_mother_mode_changes,
            "initial_charges":first_charges,"final":final,
            "max_global_charge_drift":max_initial_drift,
            "steps":rows,
            "maximum_energy_remap_defect":max(abs(r["finite_volume_remap_energy_error"]) for r in rows),
            "maximum_physical_shell_collision_energy_defect":
                max(abs(r["collision_physical_shell_energy_error"]) for r in rows),
            "time_step_convergence_not_established":True,
            "not_a_self_consistent_cosmological_model":True,
            "warning":"Prescribed RD H and constant entropy, tiny x interval, severe thermal UV truncation. Physical energy moments differ from discrete collision energies; remapping introduces an explicit error. No fully conserved FRW energy/SM bath, late freezeout, condensate or CLASS PSD."
        }


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="outputs/narrow_closed_frw_split.json")
    options=p.parse_args()
    outcome=NarrowClosedFRWPilot().run()
    path=Path(options.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(outcome,indent=2),encoding="utf-8")
    print(json.dumps(outcome,indent=2))
