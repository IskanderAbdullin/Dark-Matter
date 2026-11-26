"""Closed three-species resonance network: exact discrete invariants and H theorem."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from closed_three_species_lattice import ClosedResonantDecay


@pytest.mark.parametrize("nq,delta",[(30,.08),(50,.08),(80,.08)])
def test_two_charges_full_energy_and_entropy(nq,delta):
    k=ClosedResonantDecay(nq,delta)
    c=k.checks()
    assert c["thermal_equilibrium_residual"]<1e-11
    assert max(abs(c[x]) for x in (
        "conserved_charge1_derivative","conserved_charge2_derivative",
        "total_energy_derivative","dH_plus_sigma"))<1e-9
    assert c["entropy_production"]>0
    assert c["negative_entropy_derivative"]<0


def test_every_active_three_species_channel_has_resonant_energy_and_real_parent():
    k=ClosedResonantDecay(50,.08)
    i,j=np.nonzero(k.K)
    assert np.allclose(k.parent_E[i+j],k.e[i]+k.e[j],atol=1e-13)
    assert np.min(4*k.e[i]*k.e[j]-k.M1_T*k.M1_T)>=-1e-12
    assert np.all(k.v1[i+j]>0)


def test_all_quantum_boundaries_have_inward_derivatives():
    k=ClosedResonantDecay(30,.08)
    n=k.nq
    y=np.r_[np.full(n,.3),np.full(n,.3),np.full(2*n-1,.4)]
    z=y.copy();z[:n]=0.
    assert np.min(k.rhs(0,z)[:n])>=-1e-12
    z=y.copy();z[n:2*n]=0.
    assert np.min(k.rhs(0,z)[n:2*n])>=-1e-12
    z=y.copy();z[n:2*n]=1.
    assert np.max(k.rhs(0,z)[n:2*n])<=1e-12
    z=y.copy();z[2*n:]=0.
    assert np.min(k.rhs(0,z)[2*n:])>=-1e-12
    z=y.copy();z[2*n:]=1.
    assert np.max(k.rhs(0,z)[2*n:])<=1e-12


def test_parent_is_really_depleted_with_no_external_bath():
    k=ClosedResonantDecay(50,.08)
    result=k.evolve(tend=.03)
    initial=result["initial"]
    final=result["final"]
    assert final["n_chi1"]<initial["n_chi1"]
    assert final["n_phi"]>initial["n_phi"]
    assert final["n_chi2"]>initial["n_chi2"]
    for v in result["conservation_residuals"].values():
        assert abs(v)<1e-8
    assert result["min_occupation"]>=-1e-8
    assert result["max_fermion_occupation_chi1"]<=1.000001
    assert result["max_fermion_occupation_chi2"]<=1.000001
