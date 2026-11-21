"""Sparse, analytically differentiated CLOSED 3-species Bose/Pauli kinetics.

Extends AdaptiveClosedDecay but stores only kinematically active reaction
triplets (phi_i,chi2_j,chi1_k). A unique mother energy/occupation is shared
across every resonant pair with the same k. The collision Jacobian is
assembled analytically in CSR sparse form for SciPy's implicit BDF method.

This is a structural finite-grid solver, NOT a converged physical frozen
spectrum or a self-consistent FRW dark-matter calculation. The graph is
audited independently: default scientific mass=0.1, qmax=5 network has
three collision invariants but still 0.3% initial mother-number bias.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import coo_matrix

from adaptive_resonant_lattice import AdaptiveClosedDecay
from resonant_graph_repair import bridged_indices


class SparseClosedDecay(AdaptiveClosedDecay):
    def __init__(self, indices, delta=.001, M1_T=.1, g1=4, gchi=4):
        super().__init__(indices,delta=delta,M1_T=M1_T,g1=g1,gchi=gchi)
        self.i,self.j=np.nonzero(self.K)
        self.k=self.parent_index[self.i,self.j]
        self.w=self.K[self.i,self.j]
        self.size=2*self.nq+self.nparent
        assert np.all(self.w>0)
        assert len(self.w)==self.active_pairs
        assert self.graph["active_graph_collision_invariants"]==3,(
            "Spurious invariant graph cannot be used for this sparse study")
        assert self.graph["isolated_boson_modes"]==0
        assert self.graph["isolated_chi2_modes"]==0
        assert self.graph["inactive_declared_parent_modes"]==0

    def rhs_sparse(self,t,y):
        n=self.nq
        b,f,g=self.unpack(y)
        i,j,k=self.i,self.j,self.k
        # Exact algebraic simplification of the Bose/Pauli gain-loss:
        # g*(1-f)*(1+b) - (1-g)*f*b = g*(1+b-f) - b*f.
        r=self.w*(g[k]*(1+b[i]-f[j])-b[i]*f[j])
        return np.r_[np.bincount(i,weights=r,minlength=n)/self.vol_phi,
                     np.bincount(j,weights=r,minlength=n)/self.vol_chi,
                     -np.bincount(k,weights=r,minlength=self.nparent)/self.vol_parent]

    def jacobian(self,t,y):
        """Exact sparse Jacobian: 9 contributions per active 1<->2 triplet.

        For flux J=K[g_k(1+b_i-f_j)-b_i f_j]:
        dJ/db_i=K(g_k-f_j), dJ/df_j=-K(g_k+b_i),
        dJ/dg_k=K(1+b_i-f_j).
        """
        n=self.nq
        b,f,g=self.unpack(y)
        i,j,k=self.i,self.j,self.k
        d_b=self.w*(g[k]-f[j])
        d_f=-self.w*(g[k]+b[i])
        d_g=self.w*(1+b[i]-f[j])
        rows=[];cols=[];values=[]
        for row,scale in (
            (i,1/self.vol_phi[i]),
            (n+j,1/self.vol_chi[j]),
            (2*n+k,-1/self.vol_parent[k])
        ):
            for col,dJ in ((i,d_b),(n+j,d_f),(2*n+k,d_g)):
                rows.append(row);cols.append(col);values.append(scale*dJ)
        result=coo_matrix(
            (np.concatenate(values),
             (np.concatenate(rows),np.concatenate(cols))),
            shape=(self.size,self.size)).tocsr()
        result.sum_duplicates()
        return result

    def jacobian_directional_error(self,y,vector,epsilon=1e-6):
        y=np.asarray(y,dtype=float);vector=np.asarray(vector,dtype=float)
        if y.shape!=(self.size,) or vector.shape!=(self.size,):
            raise ValueError("Invalid Jacobian probe shapes")
        numerical=(self.rhs_sparse(0,y+epsilon*vector)-
                   self.rhs_sparse(0,y-epsilon*vector))/(2*epsilon)
        exact=self.jacobian(0,y)@vector
        return float(np.max(abs(exact-numerical))/
                     max(1.,float(np.max(abs(exact)))))

    def evolve_sparse(self,tend=.003,rtol=1e-6,atol=1e-11,
                      max_step=.0003, initial=None):
        n=self.nq
        if initial is None:
            from scipy.special import expit
            y0=np.r_[np.zeros(2*n),expit(-self.parent_E)]
        else:
            y0=np.asarray(initial,dtype=float)
        if y0.shape!=(self.size,) or tend<=0 or max_step<=0:
            raise ValueError("Invalid integration input")
        before=self.moments(y0)
        sol=solve_ivp(self.rhs_sparse,(0.,tend),y0,
                      method="BDF",jac=self.jacobian,rtol=rtol,atol=atol,
                      max_step=max_step)
        if not sol.success:
            raise RuntimeError(sol.message)
        final=sol.y[:,-1]
        b,f,g=self.unpack(final)
        after=self.moments(final)
        keys=("charge_chi1_plus_phi","charge_chi1_plus_chi2",
              "total_discrete_energy")
        absolute={key:float(after[key]-before[key]) for key in keys}
        relative={key:float(abs(absolute[key])/
                      max(1e-12,abs(before[key]))) for key in keys}
        if max(abs(v) for v in absolute.values())>1e-7:
            raise AssertionError("Conservation invariant drift")
        if min(final)<-1e-7 or max(f)>1.000001 or max(g)>1.000001:
            raise AssertionError("Quantum-statistical bounds violated")
        number_phi=float(after["n_phi"])
        return {
          "t_Gamma":float(tend),
          "nq":n,"parent_states":self.nparent,"active_reactions":self.active_pairs,
          "graph_extra_collision_invariants":self.graph[
              "additional_invariants_beyond_number_number_energy"],
          "initial":before,"final":after,
          "conserved_quantity_absolute_drift":absolute,
          "conserved_quantity_relative_drift":relative,
          "n_phi_over_T3":number_phi,
          "mean_phi_energy_over_T":float(
              np.dot(self.e*self.vol_phi,b)/number_phi) if number_phi>0 else None,
          "fraction_phi_q_below_0p1":float(
              np.dot((self.vol_phi*b)[self.e<.1],
                     np.ones(np.count_nonzero(self.e<.1)))/number_phi)
                     if number_phi>0 else None,
          "maximum_boson_occupation_grid_dependent":float(max(b)),
          "max_chi1_occupation":float(max(g)),
          "max_chi2_occupation":float(max(f)),
          "minimum_occupation":float(min(final)),
          "nfev":int(sol.nfev),"njev":int(sol.njev),"nlu":int(sol.nlu),
          "initial_mother_inventory_relative_bias":
              float(self.mother_inventory_relative_bias),
          "classification":"STRUCTURAL NUMERICAL SOLUTION; NOT frozen PSD, not CLASS ready",
        }


def run(quick=False):
    # Quick structural test at M/T=0.02: the 0.3 energy span contains
    # all kinematically required daughter IR channels, but the parent
    # thermal tail is badly truncated. Do NOT treat it as a physical
    # production rate for the main M/T=0.1 scenario.
    if quick:
        mass,maxindex,tend=.02,300,.003
    else:
        mass,maxindex,tend=.1,5000,.001
    m=SparseClosedDecay(bridged_indices(max_index=maxindex),M1_T=mass)
    rng=np.random.default_rng(913)
    y=np.r_[.1*np.exp(-m.e),.2*np.exp(-m.e),
            1/(1+np.exp(m.parent_E))]
    vec=rng.normal(size=m.size)
    jacerror=m.jacobian_directional_error(y,vec)
    rhs_error=float(np.max(abs(m.rhs(0,y)-m.rhs_sparse(0,y))))
    assert rhs_error<1e-7 and jacerror<1e-6,(rhs_error,jacerror)
    data=m.evolve_sparse(tend=tend,rtol=1e-5 if not quick else 1e-6)
    return {
      "scope":"Sparse Jacobian CLOSED 3-species 1<->2 finite-grid benchmark",
      "mass_over_T":mass,"daughter_max_index":maxindex,
      "rhs_consistency_max_abs":rhs_error,
      "jacobian_directional_relative_inf_error":jacerror,
      "evolution":data,
      "caution":"No self-consistent FRW, temperature dynamics, condensate/zero mode or true full spectral convergence; mass=.02 quick mode has severe high-energy mother truncation.",
    }


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--quick",action="store_true")
    p.add_argument("--output",default="outputs/sparse_closed_kinetics.json")
    args=p.parse_args()
    output=run(quick=args.quick)
    path=Path(args.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(output,indent=2),encoding="utf-8")
    print(json.dumps(output,indent=2))
