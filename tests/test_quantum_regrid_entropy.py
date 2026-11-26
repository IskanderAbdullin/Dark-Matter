"""Discrete quantum concavity under number-conserving daughter UV extension."""
from __future__ import annotations
import math
import sys
from pathlib import Path

import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from quantum_regrid_entropy import quantum_entropy,verify_entropy_mixing,one_full_repair


def test_bose_high_occupation_coarse_graining_has_nonnegative_entropy():
    old=[0.,.001,.01,.1]
    new=[0.,.003,.02,.1,.2]
    f=[1e8,1e3,1.]
    d=verify_entropy_mixing(old,new,f,"bose",1)
    assert d["increase"]>=-1e-8
    assert math.isfinite(d["after"])
    assert d["quantum_statistics"]=="bose"


def test_fermi_quantum_entropy_with_pauli_bounds_is_concave():
    old=[0.,.1,.3,.9]
    new=[0.,.02,.2,.5,.9,2.]
    f=[1.,0.,.4]
    d=verify_entropy_mixing(old,new,f,"fermi",4)
    assert d["increase"]>=-1e-11
    assert d["after"]>=d["before"]-1e-11


def test_pure_vacuum_has_zero_entropy_even_with_uv_extension():
    a=quantum_entropy([0,.1,1.],[0,0],"bose",1)
    b=quantum_entropy([0,.2,2.],[0,0],"fermi",4)
    assert a==0 and b==0


def test_parent_and_both_daughters_keep_two_charges_and_have_nonnegative_mixing_entropy():
    r=one_full_repair()
    assert r["two_charge_conservation_residual"]<1e-10
    assert r["total_quantum_coarse_grained_entropy_increase"]>=-1e-9
    assert all(v["increase"]>=-1e-9
               for v in r["individual_species_entropy"].values())
    assert r["numerical_entropy_not_physical_collision_production"]
    assert r["not_a_validated_cosmological_distribution"]
    assert math.isfinite(r["pure_regridding_physical_energy_error"])


def test_fermion_pauli_violation_and_invalid_statistics_are_rejected():
    with pytest.raises(ValueError):
        quantum_entropy([0,1],[1.01],"fermi",1)
    with pytest.raises(ValueError):
        quantum_entropy([0,1],[1.],"fictional",1)
