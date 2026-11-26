"""Closed multi-species number conservation when daughter UV domain must expand."""
from __future__ import annotations
import sys
from pathlib import Path
import math
import numpy as np
import pytest
from scipy.special import expit
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from daughter_uv_regrid import positive_species_regrid,number,massless_rho,remap_three_species,run
from sparse_kinetic_jacobian import SparseClosedDecay
from resonant_graph_repair import bridged_indices


def test_large_bose_occupations_are_not_spuriously_pauli_capped():
    old=[0.,.001,.01,.1]
    new=[0.,.005,.05,.1,.2]
    b=np.array([1e9,1e3,.1])
    newb=positive_species_regrid(old,new,b,"bose",1)
    assert newb.max()>1
    assert newb.min()>=0
    assert abs(number(old,b,1)-number(new,newb,1))<1e-11


def test_fermion_occupations_remain_in_unit_interval_under_uv_extension():
    old=[0.,.1,.3,1.]
    new=[0.,.05,.2,.6,1.,2.]
    f=[1.,.4,0.]
    mapped=positive_species_regrid(old,new,f,"fermi",4)
    assert 0<=min(mapped)<=max(mapped)<=1
    assert mapped[-1]==0
    assert abs(number(old,f,4)-number(new,mapped,4))<1e-12


def test_uv_domain_cannot_shrink_or_fermion_pauli_violate():
    with pytest.raises(ValueError,match="Shrinking UV"):
        positive_species_regrid([0,1,2],[0,1],[.1,.2])
    with pytest.raises(ValueError,match="Pauli"):
        positive_species_regrid([0,1],[0,2],[1.01],"fermi")
    with pytest.raises(ValueError):
        positive_species_regrid([0,1],[0,2],[-1],"bose")


def test_full_three_species_grid_change_preserves_both_comoving_charges():
    m0=SparseClosedDecay(bridged_indices(max_index=300),M1_T=.02)
    m1=SparseClosedDecay(bridged_indices(max_index=500),M1_T=.03)
    initial=np.r_[2*np.exp(-m0.e/.04),.1*np.exp(-m0.e),expit(-m0.parent_E)]
    result=remap_three_species(m0,m1,initial,.02,.03)
    final=result["remapped_state"]
    assert len(final)==2*m1.nq+m1.nparent
    assert m0.graph["active_graph_collision_invariants"]==3
    assert m1.graph["active_graph_collision_invariants"]==3
    assert result["daughter_modes_after"]>result["daughter_modes_before"]
    assert result["maximum_charge_drift"]<1e-10
    assert result["ode_moment_charge_consistency_residual"]<1e-10
    assert result["both_models_parent_phase_volumes_rebased_on_q3"]
    for name in ("phi","chi1","chi2"):
        assert abs(result["number_before"][name]-result["number_after"][name])<1e-11
    assert result["physical_free_expansion_energy_change"]>0
    assert math.isfinite(result["pure_regridding_physical_energy_error"])
    assert result["not_an_energy_conserving_remap"]


def test_cli_marks_new_high_energy_support_as_unverified_vacuum():
    result=run()
    assert result["added_momentum_domain_assumed_vacuum"]
    assert result["not_CLASS_ready"]
    assert result["remapped_occupations_length"]>400
    assert result["max_boson_occupation"]>1
    collision=result["new_graph_post_remap_collision_test"]
    assert collision["quantum_bounds_pass"]
    assert collision["nfev"]>0
    assert max(abs(v) for v in collision["discrete_charge_and_energy_residuals"].values())<1e-8
