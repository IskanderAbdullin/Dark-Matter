"""Exact 2to2 real-φ near-rest global thermal FD/BE pointwise gate."""
import math
import numpy as np
from kinetics.phi_bose_pair_nearrest_thermal_gate import (
    minimum_pair_total_energy_for_scalar_below,
    thermal_pointwise_pair_nearrest_bound,
    log10_one_plus_bose,log10_bose,
)
from kinetics.phi_bose_physical_cutcell_weak import M1,M2,MPHI
from kinetics.phi_bose_pair_ir_parent_support import (
    global_pair_scalar_support_floor,
)


def test_exact_pair_nearrest_requires_TeV_incoming_total_even_for_86p1eV():
    for m,Erest in ((M2,9302.325581395347),
                    (M1,232558.13953488387)):
        x=minimum_pair_total_energy_for_scalar_below(MPHI,m)
        assert math.isclose(x,Erest,rel_tol=5e-14)
        less=minimum_pair_total_energy_for_scalar_below(86.1e-9,m)
        assert less<x
        if m==M2:
            assert math.isclose(less,8864.413877501365,rel_tol=5e-14)
    assert minimum_pair_total_energy_for_scalar_below(86.1e-9,M2)<(
        minimum_pair_total_energy_for_scalar_below(86.1e-9,M1))


def test_exact_threshold_is_inverted_from_physical_hard_parent_pair_support():
    for m in (M2,M1):
        for K in (86.00001e-9,86.01e-9,86.1e-9,100e-9,186e-9):
            E=minimum_pair_total_energy_for_scalar_below(K,m)
            p=math.sqrt((E/2)**2-m*m)
            assert math.isclose(global_pair_scalar_support_floor(p,m),
                                K,rel_tol=1e-8)
            assert global_pair_scalar_support_floor(p*.99,m)>K
            assert global_pair_scalar_support_floor(p*1.01,m)<K


def test_subcritical_FD_BE_pair_nearrest_1GeV_pointwise_suppression():
    x=thermal_pointwise_pair_nearrest_bound(86.1e-9,M2,1.,1.)
    assert math.isclose(x["incoming_pair_total_energy_min_GeV"],
                        8864.413877501365,rel_tol=3e-13)
    assert math.isclose(x["positive_forward_log10_occupation_upper"],
                        -3842.7005307380196,rel_tol=5e-12)
    assert x["positive_reverse_log10_occupation_upper"]<-3800
    y=thermal_pointwise_pair_nearrest_bound(86.1e-9,M1,1.,1.)
    assert y["positive_forward_log10_occupation_upper"]<-95000
    z=thermal_pointwise_pair_nearrest_bound(86.1e-9,M2,100.,100.)
    assert z["positive_forward_log10_occupation_upper"]<-29


def test_thermal_high_T_caveat_and_Bose_noncritical_constraint():
    for T in (.1,1.,100.,1000.):
        r=thermal_pointwise_pair_nearrest_bound(86.1e-9,M2,T,T)
        assert math.isfinite(r["positive_forward_log10_occupation_upper"])
        assert math.isfinite(r["positive_reverse_log10_occupation_upper"])
    warm=thermal_pointwise_pair_nearrest_bound(86.1e-9,M2,1000.,1000.)
    assert warm["positive_forward_log10_occupation_upper"]>6.
    # Not a probability: stimulated Bose occupation factors can exceed 1.
    # At common chemical and thermal equilibrium F==R eventwise, so
    # this large forward bound is NOT proof of net Bose cold production.
    nearcrit=thermal_pointwise_pair_nearrest_bound(
        86.1e-9,M2,1.,1.,mu_phi=MPHI-1e-16)
    assert nearcrit["positive_forward_log10_occupation_upper"]<-3830
    with np.testing.assert_raises(ValueError):
        thermal_pointwise_pair_nearrest_bound(86.1e-9,M2,1.,1.,mu_phi=MPHI)


def test_stable_logarithms_for_tiny_subcritical_Bose_eta():
    for eta in (1e-16,1e-8,1e-4,1.,10.,100.):
        x=log10_one_plus_bose(eta)
        y=log10_bose(eta)
        assert math.isfinite(x) and math.isfinite(y)
        assert math.isclose(y,-eta/math.log(10)+x,rel_tol=1e-14,abs_tol=1e-14)
