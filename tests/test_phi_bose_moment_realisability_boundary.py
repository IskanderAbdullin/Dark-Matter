"""Finite-time affine Bose inverse-fugacity failure versus physical positive spectra.

These are frozen-source moment checks only; do not claim a full collision step.
"""
import math
import numpy as np
from kinetics.phi_bose_moment_realisability_boundary import (
    boundary_crossing, dual_gap_at_fixed_right,
    positive_constructive_update,difference_moments,
)


def test_exact_frozen_source_boundary_vs_linearized_limit():
    expected={24:174.2362815813703,
              32:26.384661164015654,
              40:149.5283248384145}
    for n,t in expected.items():
        b=boundary_crossing(n)
        assert math.isclose(b['boundary_time'],t,rel_tol=1e-7)
        assert b['two_moment_relative_crossing_residual']<1e-8
        assert .995<b['time_ratio_to_linear']<1.005
        assert b['eta_right_at_left_zero']>b['upper']


def test_affine_bose_moment_infeasible_after_boundary():
    for n in (24,32,40):
        t=boundary_crossing(n)['boundary_time']
        before=dual_gap_at_fixed_right(n,.8*t)
        after=dual_gap_at_fixed_right(n,1.1*t)
        assert before['target_minus_max_left_M']<0
        assert after['target_minus_max_left_M']>1e-14
        assert abs(after['target_right_residual'])<5e-22


def test_construct_positive_occupation_outside_affine_dual_family():
    for n in (24,32,40):
        t=boundary_crossing(n)['boundary_time']
        for k in (1.1,2.,10.):
            c=positive_constructive_update(n,k*t)
            assert c['positive']
            assert min(c['multipliers_at_edges'])>.99
            assert c['max_moment_residual']<5e-21
            assert all(v>0 for v in c['target_M'])


def test_parent_cutoff_is_essential_and_true_support_is_not_condensate():
    assert boundary_crossing(48)['support'] is False
    assert boundary_crossing(64)['support'] is False
    for pmax in (8.,16.):
        for n in (16,24,32,40):
            r=boundary_crossing(n,pmax)
            assert r['support']
            assert r['deta_left']>0
            assert r['boundary_time'] is None


def test_two_moment_boundary_exact_energy_geometry():
    b=boundary_crossing(32)
    d=difference_moments(b['upper'],0.,b['eta_right_at_left_zero'])
    C=np.array(b['C_left_right'])
    assert np.max(abs(d-b['boundary_time']*C))<5e-21
    assert b['eta_right_at_left_zero']>b['upper']
    # Finite p²dp Bose number at ηL=0 despite divergent b(p→0).
    # This is not a proof of a dynamical condensate or a physical strong step.
    assert 1.6e-7<2*86e-9*(b['upper']-86e-9)/b['eta_right_at_left_zero']<1.9e-7
