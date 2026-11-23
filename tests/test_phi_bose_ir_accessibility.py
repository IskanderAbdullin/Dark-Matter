"""Independent kinematic negative controls for 86-eV scalar IR accessibility."""
import math
import pytest
from kinetics.phi_bose_ir_accessibility import (
    M1, M2, MPHI, rest_kinematics,
    parent_threshold_for_daughter_energy, daughter_min_at_parent,
    min_over_parent_cutoff,
)


def test_exact_boosted_threshold_and_parent_rest_limit():
    estar,qstar=rest_kinematics()
    assert math.isclose(estar, .048,rel_tol=1e-10)
    assert math.isclose(parent_threshold_for_daughter_energy(estar),0,abs_tol=0)
    assert math.isclose(daughter_min_at_parent(0),estar,rel_tol=1e-14)
    pstar=parent_threshold_for_daughter_energy(MPHI)
    assert 55813.9<pstar<55814.
    assert math.isclose(daughter_min_at_parent(pstar),MPHI,rel_tol=2e-5)


@pytest.mark.parametrize("k,p_expected",[
    (MPHI,55813.9534883),
    (2*MPHI,14955.3037634),
    (10*MPHI,2797.70952099),
    (1e-6,2404.45408546),
    (1e-5,240.004427348),
    (1e-4,23.999900271),
    (1e-3,2.39895833777),
])
def test_independent_threshold_numeric_values_and_boost(k,p_expected):
    p=parent_threshold_for_daughter_energy(k)
    assert math.isclose(p,p_expected,rel_tol=2e-10)
    assert math.isclose(daughter_min_at_parent(p),k,rel_tol=2e-8)


@pytest.mark.parametrize("p,k_expected",[
    (2.,0.00119925093757947),
    (8.,0.000299988288329003),
    (16.,0.000149998547511642),
    (24.,0.0000999995844660669),
])
def test_finite_parent_support_floor(p,k_expected):
    assert math.isclose(min_over_parent_cutoff(p),k_expected,rel_tol=2e-9)
    assert min_over_parent_cutoff(p)>MPHI
    assert min_over_parent_cutoff(p)>min_over_parent_cutoff(2*p)


def test_no_below_mass_and_nonmonotone_fixed_parent_energy():
    with pytest.raises(ValueError):
        parent_threshold_for_daughter_energy(.99*MPHI)
    pstar=parent_threshold_for_daughter_energy(MPHI)
    assert min_over_parent_cutoff(2*pstar)==MPHI
    assert daughter_min_at_parent(2*pstar)>MPHI
    assert math.isclose(min_over_parent_cutoff(pstar),MPHI,rel_tol=2e-5)
