"""Frozen massive-parent distribution versus artificial thermal replenishment."""
from __future__ import annotations
import math
import sys
from pathlib import Path
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from frw_collisionless_mother import frozen_mother_background,run
from frw_reservoir_budget import thermal_mother_number


@pytest.mark.parametrize("x",[.05,.1,.5,1.,2.,5.,10.,20.])
def test_free_streaming_mother_number_and_exact_FRW_continuity(x):
    r=frozen_mother_background(x,x_dec=.05)
    assert math.isclose(r["comoving_mother_number_over_T3"],
                        thermal_mother_number(.05),rel_tol=1e-12)
    assert abs(r["FRW_collisionless_continuity_residual"])<1e-11
    assert abs(r["physical_d_ln_rho_by_d_ln_a"]-
               r["FRW_expected_d_ln_rho_by_d_ln_a"])<1e-11
    assert 0<r["equation_of_state_P_over_rho"]<1/3


def test_independent_finite_difference_of_energy_as_mass_to_T_increases():
    x=2.
    h=1e-4
    low=frozen_mother_background(x-h)
    high=frozen_mother_background(x+h)
    center=frozen_mother_background(x)
    derivative=(high["dimensionless_rho_over_T4"]-
                low["dimensionless_rho_over_T4"])/(2*h)
    assert math.isclose(derivative,
                        center["d_dimensionless_rho_dx"],
                        rel_tol=1e-9,abs_tol=1e-10)


def test_uncoupled_mother_must_not_reset_to_instantaneous_thermal_distribution():
    initial=frozen_mother_background(.05)
    late=frozen_mother_background(5.)
    later=frozen_mother_background(10.)
    assert math.isclose(initial["instantaneous_thermal_to_frozen_number_ratio"],
                        1.,rel_tol=1e-12)
    assert .07<late["instantaneous_thermal_to_frozen_number_ratio"]<.08
    assert later["instantaneous_thermal_to_frozen_number_ratio"]<.002
    assert math.isclose(initial["comoving_mother_number_over_T3"],
                        later["comoving_mother_number_over_T3"],rel_tol=1e-12)
    assert later["equation_of_state_P_over_rho"]<.04


def test_collisionless_physical_scaling_interpolates_between_radiation_and_matter():
    early=frozen_mother_background(.05)
    late=frozen_mother_background(20.)
    assert -4<early["physical_d_ln_rho_by_d_ln_a"]<-3.9
    assert -3.05<late["physical_d_ln_rho_by_d_ln_a"]<-3
    assert late["dimensionless_rho_over_T4"]>early["dimensionless_rho_over_T4"]


def test_full_background_series_is_consistent():
    result=run()
    assert len(result["cases"])==8
    assert result["cases"][0]["initial_decoupling_x"]==.05
    assert "not a" not in result["important_result"]
    assert "Does not include 1<->2 reactions" in result["warning"]


@pytest.mark.parametrize("x,initial",[(0.,.05),(-1.,.05),(1.,0.),(1.,-1.)])
def test_invalid_ratios_rejected(x,initial):
    with pytest.raises(ValueError):
        frozen_mother_background(x,x_dec=initial)
