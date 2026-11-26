"""Implicit sparse Jacobian closed quantum kinetics: dynamics and conservation."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.special import expit

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from sparse_kinetic_jacobian import SparseClosedDecay
from resonant_graph_repair import bridged_indices


@pytest.mark.parametrize("mass,maximum",[(.02,300),(.1,5000)])
def test_sparse_operator_is_identical_to_original_nonlinear_operator(mass,maximum):
    m=SparseClosedDecay(bridged_indices(max_index=maximum),M1_T=mass)
    y=np.r_[.1*np.exp(-m.e),.2*np.exp(-m.e),expit(-m.parent_E)]
    baseline=m.rhs(0,y)
    sparse=m.rhs_sparse(0,y)
    assert np.allclose(baseline,sparse,rtol=1e-12,atol=1e-8)
    n=m.nq
    dphi=np.dot(m.vol_phi,sparse[:n])
    dchi=np.dot(m.vol_chi,sparse[n:2*n])
    dparent=np.dot(m.vol_parent,sparse[2*n:])
    assert abs(dphi+dparent)<1e-10
    assert abs(dchi+dparent)<1e-10
    energy=(np.dot(m.e*m.vol_phi,sparse[:n])+
            np.dot(m.e*m.vol_chi,sparse[n:2*n])+
            np.dot(m.parent_E*m.vol_parent,sparse[2*n:]))
    assert abs(energy)<1e-10


@pytest.mark.parametrize("mass,maximum",[(.02,300),(.1,5000)])
def test_analytic_sparse_jacobian_matches_independent_finite_difference(mass,maximum):
    m=SparseClosedDecay(bridged_indices(max_index=maximum),M1_T=mass)
    rng=np.random.default_rng(2026)
    y=np.r_[.01*np.exp(-m.e),.1*np.exp(-m.e),expit(-m.parent_E)]
    v=rng.normal(size=m.size)
    error=m.jacobian_directional_error(y,v,epsilon=1e-6)
    assert error<1e-6,error
    mat=m.jacobian(0,y)
    assert mat.shape==(m.size,m.size)
    assert mat.nnz>m.size
    assert mat.nnz<m.size*m.size*.9


def test_sparse_bdf_evolution_preserves_all_conserved_quantities():
    m=SparseClosedDecay(bridged_indices(max_index=300),M1_T=.02)
    result=m.evolve_sparse(tend=.003)
    assert result["graph_extra_collision_invariants"]==0
    assert result["n_phi_over_T3"]>0
    assert result["final"]["n_chi1"]<result["initial"]["n_chi1"]
    assert result["max_chi1_occupation"]<=1.000001
    assert result["max_chi2_occupation"]<=1.000001
    assert result["minimum_occupation"]>=-1e-7
    assert result["njev"]>=1
    assert result["nlu"]>=1
    assert max(result["conserved_quantity_relative_drift"].values())<1e-7


def test_sparse_model_rejects_unconnected_unphysical_reaction_graph():
    from collision_invariant_rank import three_block_indices
    with pytest.raises(AssertionError,match="Spurious invariant graph"):
        SparseClosedDecay(three_block_indices(),M1_T=.1)
