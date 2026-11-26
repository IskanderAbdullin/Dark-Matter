"""Regression suite for the entropy-balanced finite-grid surrogate."""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from lattice_detailed_balance import EntropyBalancedLattice, scenario


@pytest.mark.parametrize("mu", [0.0,-0.3,-2.0])
def test_discrete_thermal_equilibrium(mu):
    m=EntropyBalancedLattice(50,64,48)
    y=m.equilibrium(mu)
    d=m.rhs(0,y)
    assert np.max(np.abs(d[:m.nq]*m.vol_b))<1e-11
    assert np.max(np.abs(d[m.nq:]*m.vol_f))<1e-11


def test_number_invariant_bose_positivity_fermion_pauli():
    m=EntropyBalancedLattice(50,64,48)
    x=np.r_[np.exp(-m.e),0.4*np.exp(-m.e)]
    d=m.rhs(0,x)
    assert abs(np.dot(d[:m.nq],m.vol_b)-np.dot(d[m.nq:],m.vol_f))<1e-11
    z=np.zeros(2*m.nq)
    assert np.min(m.rhs(0,z))>=-1e-12
    fone=z.copy()
    fone[m.nq:]=1.0
    assert np.max(m.rhs(0,fone)[m.nq:])<=1e-12


@pytest.mark.parametrize("factor",[0.1,1.0,10.0])
def test_discrete_entropy_production_matches_minus_free_energy(factor):
    m=EntropyBalancedLattice(50,64,48)
    x=np.r_[factor*np.exp(-m.e/0.2),0.35*np.exp(-m.e)]
    sigma=m.entropy_production(x)
    dF=m.free_energy_rate(x)
    assert sigma>0 and dF<0
    assert abs(dF+sigma)<max(1e-9,1e-10*sigma)


def test_actual_parent_quad_energy_not_conserved_by_surrogate():
    coarse=EntropyBalancedLattice(50,64,48)
    fine=EntropyBalancedLattice(100,140,72)
    assert 0<fine.absolute_parent_energy_error<coarse.absolute_parent_energy_error
    assert 0<fine.absolute_forward_weight_error<coarse.absolute_forward_weight_error
    assert 0<fine.relative_vacuum_shift<coarse.relative_vacuum_shift


def test_nonlinear_evolution_positivity_and_daughter_number():
    r=scenario(50,64,48)["evolution"]
    assert r["nphi_over_T3"]>0.1
    assert abs(r["daughter_number_difference_over_T3"])<2e-7
    assert r["min_occupation"]>=-2e-8
    assert r["max_fermion_occupation"]<=1.000001


@pytest.mark.parametrize("qmin",[1e-5,1e-6,1e-7])
def test_separate_infrared_regulator(qmin):
    m=EntropyBalancedLattice(50,64,48,qmin=qmin)
    assert m.checks()["thermal_equilibrium_max_bin_rate"]<1e-11
    assert m.absolute_parent_energy_error>0
