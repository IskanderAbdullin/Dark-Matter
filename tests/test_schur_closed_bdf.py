"""Exact Schur elimination for massive closed-parent kinetic BDF systems."""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pytest
from scipy.special import expit
from scipy.sparse import eye
from scipy.sparse.linalg import spsolve

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from schur_closed_bdf import SchurBDF,solve_closed_schur,compare_solvers_quick
from sparse_kinetic_jacobian import SparseClosedDecay
from resonant_graph_repair import bridged_indices


def test_schur_factorization_matches_independent_sparse_LU():
    model=SparseClosedDecay(bridged_indices(max_index=300),M1_T=.02)
    n=model.nq
    state=np.r_[.1*np.exp(-model.e),.2*np.exp(-model.e),
                 expit(-model.parent_E)]
    solver=SchurBDF(model.rhs_sparse,0.,state,1.,n_daughters=2*n,
                    jac=model.jacobian)
    J=model.jacobian(0.,state)
    matrix=eye(model.size,format="csc")-.0001*J
    random=np.random.default_rng(271828).normal(size=model.size)
    factors=solver._factor_schur(matrix)
    actual=solver._solve_schur(factors,random)
    reference=spsolve(matrix,random)
    assert np.allclose(actual,reference,rtol=1e-10,atol=1e-9)
    assert np.max(abs(matrix@actual-random))<1e-8


def test_schur_bdf_is_same_nonlinear_equation_as_standard_implicit_BDF():
    result=compare_solvers_quick()
    assert result["absolute_n_phi_difference"]<1e-10
    a=result["ordinary_BDF"]
    b=result["Schur_BDF"]
    assert abs(a["final"]["n_chi1"]-b["n_mother_over_T3"])<1e-10
    assert b["graph_collision_invariants"]==3
    assert max(abs(r) for r in b["number_and_energy_residuals"].values())<1e-8


def test_large_physical_mass_chi1_is_dynamically_depleted_on_clean_graph():
    model=SparseClosedDecay(bridged_indices(max_index=5000),M1_T=.1)
    out=solve_closed_schur(model,tend=.0002)
    assert out["state_variables"]>5000
    assert out["parent_modes"]>5000
    assert out["n_phi_over_T3"]>0
    assert out["n_mother_over_T3"]<model.mother_thermal_number_discretized
    assert out["graph_collision_invariants"]==3
    assert out["minimum_occupation"]>=-1e-7
    assert out["max_fermion_occupation"]<=1.000001
    assert max(abs(x) for x in out["number_and_energy_residuals"].values())<1e-8


def test_invalid_schur_layout_rejected():
    model=SparseClosedDecay(bridged_indices(max_index=300),M1_T=.02)
    state=np.r_[np.zeros(2*model.nq),expit(-model.parent_E)]
    with pytest.raises(ValueError):
        SchurBDF(model.rhs_sparse,0.,state,1.,n_daughters=0,
                 jac=model.jacobian)
