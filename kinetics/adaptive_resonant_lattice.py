"""Experimental ADAPTIVE closed three-species resonant energy lattice.

Use a sparse, strictly increasing SUBSET of integer daughter energy nodes
e_i=(index_i+.5)*delta. Unique mother energies are exact pairwise sums:
E_k=e_i+e_j. A SINGLE mother occupation is shared by all pairs with
the same energy, avoiding independent replicated mother variables.
All three distributions evolve with exact two charges and discrete energy.

Crucially, a nonuniform sampled reaction graph can have SPURIOUS collision
invariants. We audit these explicitly. Numerical results from meshes that
fail the graph-rank gate MUST NOT be used as physical PSD, CLASS input,
a claimed BEC result or converged cosmological production.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expit

from closed_three_species_lattice import allowed_area_from_zero
from collision_invariant_rank import (
    audit_integer_grid, sparse_cold_hot_indices, three_block_indices
)
from frw_reservoir_budget import thermal_mother_number


class AdaptiveClosedDecay:
    def __init__(self, indices, delta=.001, M1_T=.1, g1=4, gchi=4):
        index=np.asarray(sorted(set(int(i) for i in indices)),dtype=np.int64)
        if len(index)<3 or np.any(index<0) or delta<=0 or M1_T<=0:
            raise ValueError("Invalid adaptive resonant grid")
        self.index=index;self.nq=len(index);self.delta=delta
        self.M1_T=M1_T
        self.e=(index+.5)*delta
        self.edges=np.r_[0.,(self.e[:-1]+self.e[1:])/2,
                          self.e[-1]+(self.e[-1]-self.e[-2])/2]
        self.vol_phi=(self.edges[1:]**3-self.edges[:-1]**3)/(6*math.pi**2)
        self.vol_chi=gchi*self.vol_phi
        self.graph=audit_integer_grid(index,delta=delta,mass=M1_T)
        sums=index[:,None]+index[None,:]+1
        active=4*self.e[:,None]*self.e[None,:]>=M1_T*M1_T
        parent_keys=np.unique(sums[active])
        self.parent_E=parent_keys*delta
        nparent=len(parent_keys)
        if nparent<2:raise ValueError("Too few active parent modes")
        self.nparent=nparent
        p_lo=max(M1_T,self.parent_E[0]-(self.parent_E[1]-self.parent_E[0])/2)
        parent_edges=np.r_[p_lo,(self.parent_E[:-1]+self.parent_E[1:])/2,
                            self.parent_E[-1]+(self.parent_E[-1]-self.parent_E[-2])/2]
        mom_hi=np.maximum(0.,parent_edges[1:]**2-M1_T*M1_T)**1.5
        mom_lo=np.maximum(0.,parent_edges[:-1]**2-M1_T*M1_T)**1.5
        self.vol_parent=g1/(6*math.pi**2)*(mom_hi-mom_lo)
        if np.any(self.vol_parent<=0):
            raise ValueError("Invalid mother shell volume")
        # All active pairs map to a UNIQUE shared mother occupation.
        parent_index=np.searchsorted(parent_keys,sums)
        parent_index[~active]=0
        self.parent_index=parent_index
        self.K=(g1*M1_T/(2*math.pi**2))*allowed_area_from_zero(
            self.edges,M1_T)*active
        self.active_pairs=int(np.count_nonzero(self.K))
        chosen=np.nonzero(self.K)
        self.max_resonance_error=float(np.max(abs(
            self.parent_E[self.parent_index[chosen]]-
            self.e[chosen[0]]-self.e[chosen[1]])))
        assert self.max_resonance_error<1e-12
        self.mother_thermal_number_exact=thermal_mother_number(M1_T,g1=g1)
        self.mother_thermal_number_discretized=float(np.dot(
            self.vol_parent,expit(-self.parent_E)))
        self.mother_inventory_relative_bias=(
            self.mother_thermal_number_discretized/
            self.mother_thermal_number_exact-1)
        self.not_ready_for_cosmology=True

    def unpack(self,y):
        n=self.nq
        return y[:n],y[n:2*n],y[2*n:]

    def rhs(self,t,y):
        n=self.nq
        b,f,g=self.unpack(y)
        parent=g[self.parent_index]
        J=self.K*(parent*(1-f[None,:])*(1+b[:,None])-
                  (1-parent)*f[None,:]*b[:,None])
        dparent=-np.bincount(self.parent_index.ravel(),weights=J.ravel(),
                            minlength=self.nparent)/self.vol_parent
        return np.r_[J.sum(axis=1)/self.vol_phi,
                     J.sum(axis=0)/self.vol_chi,dparent]

    def equilibrium(self,mu_phi=-.3,mu_chi=.1):
        if mu_phi>0:raise ValueError("Boson chemical potential > 0")
        return np.r_[1/np.expm1(self.e-mu_phi),
                     expit(mu_chi-self.e),
                     expit(mu_phi+mu_chi-self.parent_E)]

    def moments(self,y):
        b,f,g=self.unpack(y)
        nb=float(np.dot(self.vol_phi,b))
        nf=float(np.dot(self.vol_chi,f))
        ng=float(np.dot(self.vol_parent,g))
        en=float(np.dot(self.vol_phi*self.e,b)+
                 np.dot(self.vol_chi*self.e,f)+
                 np.dot(self.vol_parent*self.parent_E,g))
        return {"n_phi":nb,"n_chi2":nf,"n_chi1":ng,
                "charge_chi1_plus_phi":ng+nb,
                "charge_chi1_plus_chi2":ng+nf,
                "total_discrete_energy":en}

    def thermodynamic_checks(self):
        n=self.nq
        yeq=self.equilibrium()
        deq=self.rhs(0.,yeq)
        eq=float(np.max(np.abs(deq*np.r_[
            self.vol_phi,self.vol_chi,self.vol_parent])))
        y=np.r_[0.1*np.exp(-self.e),0.2*np.exp(-self.e),
                expit(-self.parent_E)]
        d=self.rhs(0.,y)
        charges=float(max(
            abs(np.dot(self.vol_phi,d[:n])+
                np.dot(self.vol_parent,d[2*n:])),
            abs(np.dot(self.vol_chi,d[n:2*n])+
                np.dot(self.vol_parent,d[2*n:]))))
        energy=float(np.dot(self.e*self.vol_phi,d[:n])+
                     np.dot(self.e*self.vol_chi,d[n:2*n])+
                     np.dot(self.parent_E*self.vol_parent,d[2*n:]))
        b,f,g=self.unpack(y)
        H=float(np.dot(self.vol_phi*d[:n],np.log(b/(1+b)))+
                np.dot(self.vol_chi*d[n:2*n],np.log(f/(1-f)))+
                np.dot(self.vol_parent*d[2*n:],np.log(g/(1-g))))
        a=g[self.parent_index]*(1-f[None,:])*(1+b[:,None])
        z=(1-g[self.parent_index])*f[None,:]*b[:,None]
        active=self.K>0
        sigma=float(np.sum((self.K*(a-z))[active]*
                           np.log(a[active]/z[active])))
        out={"equilibrium_residual":eq,
             "number_charge_derivative_max_abs":charges,
             "energy_derivative_abs":abs(energy),
             "negative_entropy_rate":H,
             "entropy_production":sigma,
             "entropy_balance_residual":H+sigma}
        assert eq<1e-11 and charges<1e-10 and abs(energy)<1e-10,out
        assert sigma>0 and H<0 and abs(H+sigma)<1e-8,out
        return out

    def evolve(self,tend=.03,rtol=1e-6,atol=1e-12):
        n=self.nq
        initial=np.r_[np.zeros(2*n),expit(-self.parent_E)]
        before=self.moments(initial)
        sol=solve_ivp(self.rhs,(0.,tend),initial,method="BDF",
                      rtol=rtol,atol=atol,max_step=.0003)
        if not sol.success:raise RuntimeError(sol.message)
        y=sol.y[:,-1];b,f,g=self.unpack(y)
        after=self.moments(y)
        invariant_keys=("charge_chi1_plus_phi","charge_chi1_plus_chi2",
                        "total_discrete_energy")
        residuals={key:after[key]-before[key] for key in invariant_keys}
        if max(abs(x) for x in residuals.values())>1e-8:
            raise AssertionError("Full energy/number budget failed")
        if min(y)<-1e-8 or max(f)>1.000001 or max(g)>1.000001:
            raise AssertionError("Quantum bounds failed")
        return {
            "t_Gamma":tend,"initial":before,"final":after,
            "charges_and_energy_residuals":residuals,
            "mean_phi_q":float(np.dot(b*self.vol_phi,self.e)/after["n_phi"]),
            "fraction_phi_q_below_0p1":float(
                np.sum((b*self.vol_phi)[self.e<.1])/after["n_phi"]),
            "min_occupation":float(min(y)),
            "max_bose_occupation_grid_dependent":float(max(b)),
            "max_mother_fermion_occupation":float(max(g)),
            "max_daughter_fermion_occupation":float(max(f)),
            "nfev":int(sol.nfev),
        }


def one_case(name,indices):
    model=AdaptiveClosedDecay(indices)
    return {
      "name":name,
      "daughter_grid_size":model.nq,
      "unique_parent_energy_modes":model.nparent,
      "active_reactions":model.active_pairs,
      "spurious_collision_invariants":model.graph[
          "additional_invariants_beyond_number_number_energy"],
      "parent_inventory_relative_bias":model.mother_inventory_relative_bias,
      "parent_initial_n_over_T3":model.mother_thermal_number_discretized,
      "exact_thermal_parent_n_over_T3":model.mother_thermal_number_exact,
      "max_resonance_energy_error":model.max_resonance_error,
      "structural_checks":model.thermodynamic_checks(),
      "evolution":model.evolve(),
      "rejected_as_physical_dark_matter_spectrum":True,
      "reason":"The graph has nonphysical collision invariants and/or mother-shell discretization bias. Numerical moment conservation alone does not establish convergence."
    }


def run(quick=False):
    cases=[one_case("sparse_cold_hot_27",sparse_cold_hot_indices())]
    if not quick:
        cases.append(one_case("three_blocks_56",three_block_indices()))
    return {"title":"Adaptive resonant CLOSED 1<->2 proof-of-concept",
            "warning":"This is a deliberately nonconverged structural experiment, NOT a primordial f(q), Omega_DM prediction or reliable BEC model.",
            "method":"Unique mother energy for every exactly resonant sum; no duplicate mother variables; analytic on-shell area weights.",
            "scenarios":cases}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--quick",action="store_true")
    parser.add_argument("--output",default="outputs/adaptive_closed_resonant.json")
    args=parser.parse_args()
    out=run(args.quick)
    path=Path(args.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,indent=2),encoding="utf-8")
    print(json.dumps(out,indent=2))
