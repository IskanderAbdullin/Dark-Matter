"""Symmetric (Strang-type) closed 1<->2 plus Liouville split in a narrow FRW test.

At fixed prescribed constant-gstar radiation background, repeat
 half collision(x_left) -> q^3 parent remap -> half collision(x_right).
Each collision is a full quantum gain-loss ODE with current exact
on-shell energy triples and a unique shared mother occupation.
The remap conserves Pauli+number, but not physical energy or entropy.
An x-dependent reaction threshold makes second-order convergence
NONTRIVIAL; compare with left-endpoint Lie splitting instead of
assuming a formal Strang order.

This is an intentionally short structurally closed toy, not a
converged physical freezeout or realistic CLASS primordial spectrum.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import numpy as np

from narrow_closed_frw_split import NarrowClosedFRWPilot
from parent_q_remap import thermal_shell_averages, rebin


class SymmetricFRWPilot(NarrowClosedFRWPilot):
    def run(self):
        m0=self.models[0]
        g=thermal_shell_averages(self.parent_edges[0],float(self.xgrid[0]))
        state=np.r_[np.zeros(2*self.n),g]
        q0=self._charges(m0,state)
        rows=[]
        for index in range(len(self.models)-1):
            left=self.models[index];right=self.models[index+1]
            x=float(self.xgrid[index]);y=float(self.xgrid[index+1])
            before=self._charges(left,state)
            half_tau=self.kappa*(y*y-x*x)/4
            discrete0=left.moments(state)["total_discrete_energy"]
            phys0=self.physical_rho(left,self.parent_edges[index],state,x)
            halfstate,ca=self._collide(left,state,half_tau)
            collision_left_discrete=left.moments(halfstate)["total_discrete_energy"]-discrete0
            collision_left_physical=self.physical_rho(
                left,self.parent_edges[index],halfstate,x)-phys0
            at_y_unchanged=self.physical_rho(left,self.parent_edges[index],halfstate,y)
            free_energy_change=at_y_unchanged-self.physical_rho(
                left,self.parent_edges[index],halfstate,x)
            b,f,g=left.unpack(halfstate)
            mapped=rebin(self.parent_edges[index],self.parent_edges[index+1],g)
            mid=np.r_[b.copy(),f.copy(),mapped]
            remap_error=self.physical_rho(
                right,self.parent_edges[index+1],mid,y)-at_y_unchanged
            discrete_mid=right.moments(mid)["total_discrete_energy"]
            phys_mid=self.physical_rho(right,self.parent_edges[index+1],mid,y)
            afterstate,cb=self._collide(right,mid,half_tau)
            collision_right_discrete=(
                right.moments(afterstate)["total_discrete_energy"]-discrete_mid)
            collision_right_physical=(
                self.physical_rho(right,self.parent_edges[index+1],afterstate,y)-
                phys_mid)
            state=afterstate
            after=self._charges(right,state)
            charge_error=max(abs(before[0]-after[0]),abs(before[1]-after[1]))
            if charge_error>1e-8:
                raise ArithmeticError("Comoving number charges not conserved")
            if max(abs(collision_left_discrete),abs(collision_right_discrete))>1e-8:
                raise ArithmeticError("Collision discrete energies not conserved")
            rows.append({
                "x_left":x,"x_right":y,"delta_tGamma":2*half_tau,
                "left_on_shell_pairs":left.active_pairs,
                "right_on_shell_pairs":right.active_pairs,
                "parent_energy_modes_before":left.nparent,
                "parent_energy_modes_after":right.nparent,
                "maximum_charge_error":charge_error,
                "left_collision_discrete_energy_error":collision_left_discrete,
                "right_collision_discrete_energy_error":collision_right_discrete,
                "left_collision_physical_shell_error":collision_left_physical,
                "right_collision_physical_shell_error":collision_right_physical,
                "physical_free_expansion_energy_change":free_energy_change,
                "numerical_remap_energy_error":remap_error,
                "left_nfev":ca["nfev"],"right_nfev":cb["nfev"]
            })
        final=self.models[-1].moments(state)
        global_error=max(
            abs(final["charge_chi1_plus_phi"]-q0[0]),
            abs(final["charge_chi1_plus_chi2"]-q0[1]))
        if global_error>1e-8:
            raise ArithmeticError("Global two-charge conservation failed")
        return {
            "method":"symmetric left/right endpoint collision halves with middle conservative q^3 parent remap",
            "model":"closed Bose-Pauli 1<->2, prescribed RD H and constant gstarS",
            "daughter_bins":self.n,"parent_modes_initial":self.models[0].nparent,
            "parent_modes_final":self.models[-1].nparent,
            "initial_x":float(self.xgrid[0]),"final_x":float(self.xgrid[-1]),
            "time_steps":len(rows),"kappa":self.kappa,
            "initial_charges":q0,"final":final,
            "global_two_charge_drift":global_error,
            "steps":rows,
            "physical_total_energy_conservation_claimed":False,
            "time_refinement_not_assumed":True,
            "CLASS_ready":False,
            "warning":"Physical shell energy bias plus nonzero numerical remap drift. Tiny x interval and truncated parent thermal tail; not full entropy/energy-closed FRW or frozen DM PSD."
        }


def compare_lie_and_symmetric(xend=.0202,kappa=500.):
    from narrow_closed_frw_split import NarrowClosedFRWPilot
    out={}
    for name,cls in (("left_endpoint_Lie",NarrowClosedFRWPilot),
                     ("symmetric",SymmetricFRWPilot)):
        cases=[]
        for steps in (2,4):
            result=cls(x0=.02,xend=xend,steps=steps,kappa=kappa).run()
            cases.append({"steps":steps,"n_phi_over_T3":result["final"]["n_phi"],
                          "charge_drift":result.get("global_two_charge_drift",
                                                    result.get("max_global_charge_drift"))})
        out[name]={"runs":cases,
                   "relative_step_refinement_difference":
                       abs(cases[0]["n_phi_over_T3"]/cases[1]["n_phi_over_T3"]-1)}
    out["symmetric_has_smaller_observed_step_error"]=bool(
        out["symmetric"]["relative_step_refinement_difference"]<
        out["left_endpoint_Lie"]["relative_step_refinement_difference"])
    out["no_physical_time_convergence_claim"]=True
    return out


def run():
    comparison=compare_lie_and_symmetric()
    modes=SymmetricFRWPilot(x0=.02,xend=.022,steps=4,kappa=10.,
                             allow_mother_mode_changes=True).run()
    return {"comparison":comparison,
            "varying_mother_energy_basis": {
                "parent_modes_initial":modes["parent_modes_initial"],
                "parent_modes_final":modes["parent_modes_final"],
                "phi_number":modes["final"]["n_phi"],
                "charge_drift":modes["global_two_charge_drift"]},
            "scientific_status":"EXPERIMENTAL, NON-COSMOLOGICAL, SHORT INTERVAL"}


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="outputs/strang_frw_comparison.json")
    opts=p.parse_args()
    out=run()
    path=Path(opts.output)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,indent=2),encoding="utf-8")
    print(json.dumps(out,indent=2))
