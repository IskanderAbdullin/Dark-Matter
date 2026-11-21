"""Adaptive resonant CLOSED three-species solver: preserve charges, reject spurious grids."""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from adaptive_resonant_lattice import AdaptiveClosedDecay
from collision_invariant_rank import sparse_cold_hot_indices,three_block_indices


@pytest.mark.parametrize("indices",[sparse_cold_hot_indices(),three_block_indices()])
def test_adaptive_unique_mother_energies_and_on_shell_reactions(indices):
    m=AdaptiveClosedDecay(indices)
    n=m.nq
    a,b=np.nonzero(m.K)
    assert m.max_resonance_error<1e-12
    assert np.allclose(m.parent_E[m.parent_index[a,b]],m.e[a]+m.e[b],
                       rtol=0,atol=1e-12)
    assert np.min(4*m.e[a]*m.e[b]-m.M1_T**2)>-1e-12
    assert np.all(m.vol_parent>0)
    assert m.nparent==m.graph["active_unique_mother_energy_modes"]
    assert m.not_ready_for_cosmology
    assert m.graph["additional_invariants_beyond_number_number_energy"]>0


@pytest.mark.parametrize("indices",[sparse_cold_hot_indices(),three_block_indices()])
def test_thermodynamic_operators_conserve_all_three_physically_required_moments(indices):
    m=AdaptiveClosedDecay(indices)
    r=m.thermodynamic_checks()
    assert r["equilibrium_residual"]<1e-11
    assert r["number_charge_derivative_max_abs"]<1e-10
    assert r["energy_derivative_abs"]<1e-10
    assert r["negative_entropy_rate"]<0
    assert r["entropy_production"]>0
    assert abs(r["entropy_balance_residual"])<1e-8


def test_adaptive_mesh_bias_is_independent_failure_mode():
    a=AdaptiveClosedDecay(sparse_cold_hot_indices())
    b=AdaptiveClosedDecay(three_block_indices())
    assert abs(a.mother_inventory_relative_bias)<.01
    assert b.mother_inventory_relative_bias<-.1
    assert a.graph["additional_invariants_beyond_number_number_energy"]>1
    assert b.graph["additional_invariants_beyond_number_number_energy"]==1


def test_cold_state_evolution_depletes_one_shared_mother_reservoir():
    model=AdaptiveClosedDecay(sparse_cold_hot_indices())
    out=model.evolve(tend=.03)
    before,after=out["initial"],out["final"]
    assert after["n_phi"]>0
    assert after["n_chi2"]>0
    assert after["n_chi1"]<before["n_chi1"]
    assert after["n_phi"]<before["n_chi1"]
    assert abs(after["n_phi"]-(before["n_chi1"]-after["n_chi1"]))<1e-8
    assert max(abs(x) for x in out["charges_and_energy_residuals"].values())<1e-8
    assert out["min_occupation"]>=-1e-8
    assert out["max_mother_fermion_occupation"]<=1.000001
    assert out["max_daughter_fermion_occupation"]<=1.000001
    assert out["max_bose_occupation_grid_dependent"]>1e6


def test_adaptive_solution_is_explicitly_rejected_as_cosmology():
    m=AdaptiveClosedDecay(three_block_indices())
    assert m.not_ready_for_cosmology
    assert not m.graph["accept_for_unique_chemical_equilibrium"]


def test_extra_invariant_is_conserved_by_full_nonlinear_adaptive_ode():
    from collision_invariant_rank import extra_invariant_coefficients
    for indices in (sparse_cold_hot_indices(),three_block_indices()):
        m=AdaptiveClosedDecay(indices)
        w=extra_invariant_coefficients(indices,delta=m.delta,mass=m.M1_T)
        assert np.array_equal(
            w["parent_integer_energy_keys"],
            np.unique((m.index[:,None]+m.index[None,:]+1)[m.K>0]))
        n=m.nq
        y=np.r_[0.1*np.exp(-m.e),0.2*np.exp(-m.e),
                1/(1+np.exp(m.parent_E))]
        dot=m.rhs(0,y)
        fake_charge_rate=(
            np.dot(w["boson_coefficients"]*m.vol_phi,dot[:n])+
            np.dot(w["chi2_coefficients"]*m.vol_chi,dot[n:2*n])+
            np.dot(w["parent_coefficients"]*m.vol_parent,dot[2*n:])
        )
        assert abs(fake_charge_rate)<1e-10
        assert w["max_reaction_residual"]<1e-8
