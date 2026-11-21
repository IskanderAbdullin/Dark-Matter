"""Independent momentum-differential source matches exact 1D thermal vacuum rate."""
from __future__ import annotations
import math
import sys
from pathlib import Path
import numpy as np
import pytest
from scipy.integrate import quad

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from best_vacuum_spectrum import (
    vacuum_phi_source_differential,analytical_peak,
    compare_best_phi_grid,run
)
from onshell_collision_2d import vacuum_reference


def test_integrating_differential_source_recovers_exact_original_number_rate():
    total=quad(lambda q:vacuum_phi_source_differential(q),0,np.inf,
               epsabs=1e-12,epsrel=1e-11)[0]
    assert math.isclose(total,vacuum_reference(),rel_tol=1e-10)
    assert math.isclose(total,.0164959084596,rel_tol=1e-10)


def test_vacuum_source_exact_peak_at_half_m_over_T():
    for m in (.05,.1,.3):
        peak=analytical_peak(m)
        assert peak["q_peak"]==m/2
        assert vacuum_phi_source_differential(m/2,m)==peak["source_max_dR_dq"]
        assert vacuum_phi_source_differential(m/4,m)<peak["source_max_dR_dq"]
        assert vacuum_phi_source_differential(m,m)<peak["source_max_dR_dq"]
    assert vacuum_phi_source_differential(0.)==0.


def test_source_is_exponentially_suppressed_in_initial_deep_ir():
    low=vacuum_phi_source_differential(1e-4)
    at_peak=vacuum_phi_source_differential(.05)
    assert low/at_peak<1e-10
    assert vacuum_phi_source_differential(1e-6)<low
    assert at_peak>0


def test_exact_synthetic_per_mode_occupation_rate_matches_spectral_reference():
    q=np.geomspace(1e-3,20,200)
    expected=vacuum_phi_source_differential(q)
    collision=2*np.pi**2*expected/q**2
    report=compare_best_phi_grid(q,collision)
    assert report["robust_points_count"]>10
    assert report["relative_peak_normalized_rms_difference"]<1e-14
    assert report["maximum_absolute_relative_difference_within_support"]<1e-14
    assert report["not_a_nonlinear_freezeout_shape_test"]


def test_wrong_occupation_rate_normalization_is_detected():
    q=np.geomspace(1e-3,10,100)
    reference=vacuum_phi_source_differential(q)
    wrong=2*np.pi**2*reference/q**2*.7
    report=compare_best_phi_grid(q,wrong)
    assert .29<report["maximum_absolute_relative_difference_within_support"]<.31


@pytest.mark.parametrize("badq",[
    [0,1,2,3,4,5,6,7,8,9],
    [1,2,3,4,5,5,7,8,9,10]])
def test_invalid_external_grid_rejected(badq):
    with pytest.raises(ValueError):
        compare_best_phi_grid(badq,np.ones(10))


def test_source_samples_are_labeled_vacuum_only_and_external_unrun():
    result=run()
    assert result["independent_BEST_executed"] is False
    assert result["peak"]["q_peak"]==.05
    assert abs(result["number_rate_integrated_from_shape"]-
               result["independent_1D_full_number_rate"])<1e-10
    assert len(result["source_samples"])>=10
