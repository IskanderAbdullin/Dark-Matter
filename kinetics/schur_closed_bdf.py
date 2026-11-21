"""Schur-complement linear solver for huge-mother-state 1<->2 BDF kinetics.

The closed kinetic Jacobian has exact block structure:
       [ Jdd  Jdp ]
  J =  [ Jpd  Jpp ]
where d=2*N_daughter, p=N_unique_parent and Jpp is diagonal.
BDF Newton solves (I-cJ)z=r. Eliminating diagonal parent block
reduces sparse factorization of ~11k variables to dense LU of only
~600 daughter variables, with no physical approximation.

The implementation subclasses SciPy's *private* BDF class to replace
the LU kernel only. Regression tests pin parity against standard BDF;
the private API can change, so this must fail closed on incompatible
SciPy versions. No modification of the physical collision operator.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate._ivp.bdf import BDF
from scipy.linalg import lu_factor,lu_solve
from scipy.sparse import issparse

from sparse_kinetic_jacobian import SparseClosedDecay
from resonant_graph_repair import bridged_indices


class SchurBDF(BDF):
    def __init__(self,fun,t0,y0,t_bound,*,n_daughters,**kwargs):
        self.n_daughters=int(n_daughters)
        if self.n_daughters<=0 or self.n_daughters>=len(y0):
            raise ValueError("must supply 0<n_daughters<len(y0)")
        super().__init__(fun,t0,y0,t_bound,**kwargs)
        if not issparse(self.J):
            raise TypeError("SchurBDF requires a sparse analytic Jacobian")
        if not hasattr(self,"lu") or not hasattr(self,"solve_lu"):
            raise RuntimeError("SciPy BDF internal LU interface changed")
        self.lu=self._factor_schur
        self.solve_lu=self._solve_schur

    def _factor_schur(self,A):
        """Factor the exact Schur complement of the diagonal parent block."""
        if not issparse(A):
            raise TypeError("Expected sparse I-cJ matrix from SciPy BDF")
        d=self.n_daughters
        P=A[d:,d:].tocoo()
        if np.any(P.row!=P.col):
            raise RuntimeError("Parent-parent Jacobian is not diagonal")
        diag=A.diagonal()[d:]
        if np.any(abs(diag)<1e-13):
            raise np.linalg.LinAlgError("Near-singular parent diagonal")
        pinv=1./diag
        B=A[:d,d:].tocsr()
        C=A[d:,:d].tocsr()
        # Exact block Gaussian elimination; sparse off-diagonal products,
        # small dense daughter Schur LU (not huge mother sparse LU).
        S=A[:d,:d].toarray()-(B.multiply(pinv)@C).toarray()
        factor=lu_factor(S,check_finite=False)
        self.nlu+=1
        return factor,B,C,pinv

    def _solve_schur(self,factors,rhs):
        factor,B,C,pinv=factors
        d=self.n_daughters
        y=rhs[:d]-B@(pinv*rhs[d:])
        daughter=lu_solve(factor,y,check_finite=False)
        parent=pinv*(rhs[d:]-C@daughter)
        return np.r_[daughter,parent]


def solve_closed_schur(model,tend=.001,rtol=1e-5,atol=1e-11,
                       max_step=.0003,return_spectrum=False):
    from scipy.special import expit
    n=model.nq
    initial=np.r_[np.zeros(2*n),expit(-model.parent_E)]
    before=model.moments(initial)
    sol=solve_ivp(model.rhs_sparse,(0.,tend),initial,method=SchurBDF,
                  n_daughters=2*n,jac=model.jacobian,
                  rtol=rtol,atol=atol,max_step=max_step)
    if not sol.success:
        raise RuntimeError(sol.message)
    y=sol.y[:,-1]
    b,f,g=model.unpack(y)
    final=model.moments(y)
    keys=("charge_chi1_plus_phi","charge_chi1_plus_chi2",
          "total_discrete_energy")
    residual={k:float(final[k]-before[k]) for k in keys}
    assert max(abs(x) for x in residual.values())<1e-7,residual
    assert min(y)>=-1e-7 and max(f)<=1.000001 and max(g)<=1.000001
    phi=final["n_phi"]
    report={
       "solver":"BDF with mathematically exact daughter Schur / diagonal mother elimination",
       "t_Gamma":tend,"daughter_modes":n,"parent_modes":model.nparent,
       "state_variables":model.size,"active_pairs":model.active_pairs,
       "nfev":int(sol.nfev),"njev":int(sol.njev),"nlu":int(sol.nlu),
       "n_phi_over_T3":float(phi),
       "n_chi2_over_T3":float(final["n_chi2"]),
       "n_mother_over_T3":float(final["n_chi1"]),
       "mean_phi_q":float(np.dot(b*model.vol_phi,model.e)/phi),
       "phi_fraction_below_q0p1":float((b*model.vol_phi)[model.e<.1].sum()/phi),
       "max_boson_occupation_grid_dependent":float(max(b)),
       "max_fermion_occupation":float(max(max(f),max(g))),
       "min_occupation":float(min(y)),
       "number_and_energy_residuals":residual,
       "relative_parent_thermal_number_error_at_initial":
           float(model.mother_inventory_relative_bias),
       "graph_collision_invariants":int(model.graph["active_graph_collision_invariants"]),
       "NOT_a_cosmological_freezeout":True,
    }
    if return_spectrum:
        report["daughter_spectrum"]={
            "q_bin_edges":list(map(float,model.edges)),
            "q_bin_representatives":list(map(float,model.e)),
            "phi_occupations":list(map(float,b)),
            "chi2_occupations":list(map(float,f)),
            "phi_number_per_bin":list(map(float,b*model.vol_phi)),
            "total_phi_number":float(phi),
            "occupation_is_grid_dependent":True,
            "spectral_status":"SHORT-TIME, NOT FROZEN, NOT VALIDATED FOR CLASS",
        }
    return report


def compare_solvers_quick():
    model=SparseClosedDecay(bridged_indices(max_index=300),M1_T=.02)
    expected=model.evolve_sparse(tend=.003,rtol=1e-5,atol=1e-11)
    schur=solve_closed_schur(model,tend=.003,rtol=1e-5,atol=1e-11)
    absolute=abs(expected["n_phi_over_T3"]-schur["n_phi_over_T3"])
    assert absolute<1e-10,(expected,schur)
    assert abs(schur["nfev"]-expected["nfev"])<=4,(expected,schur)
    return {"ordinary_BDF":expected,"Schur_BDF":schur,
            "absolute_n_phi_difference":absolute}


def run(quick=False):
    if quick:
        return {"status":"structural solver parity, truncated parent source",
                "comparison":compare_solvers_quick()}
    rows=[]
    for max_index in (5000,7000,8500,10000):
        model=SparseClosedDecay(bridged_indices(max_index=max_index))
        result=solve_closed_schur(model)
        rows.append({"max_daughter_index":max_index,"result":result})
    return {
      "status":"SHORT CLOSED 1<->2 KINETIC CUTOFF CONVERGENCE; NOT complete freeze-out",
      "source":"M1/T=.1, Delta=.001, fine core 0..99, coarse stride 50",
      "time":"tGamma=.001, fixed T flat space, exact on-shell finite-cell rates",
      "strong_conclusion":"Parent number/energy/pressure and short-time Bose production converge separately; their errors are NOT equivalent.",
      "warnings":"This is an adaptive-cell approximation with a still-unresolved zero mode and no FRW. Converged integrated number does not prove full f(q) spectral convergence or physical DM.",
      "scenarios":rows,
    }


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--quick",action="store_true")
    p.add_argument("--output",default="outputs/closed_schur_convergence.json")
    args=p.parse_args()
    res=run(quick=args.quick)
    dest=Path(args.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(res,indent=2),encoding="utf-8")
    print(json.dumps(res,indent=2))
