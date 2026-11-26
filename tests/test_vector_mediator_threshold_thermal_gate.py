"""Conditional vector-mediator equilibrium and exact on-shell threshold gates."""
import math
import pytest
from kinetics.vector_mediator_threshold_thermal_gate import (
    mediator_ideal_equilibrium,independent_bessel_eos,gauge_thresholds,dirac_vector_width_coefficient)

@pytest.mark.parametrize("mass,T",[(.03,1),(.1,1),(.3,1),(.3,.1),(3,1)])
def test_ideal_spin_one_eos_vs_independent_bessel(mass,T):
    row=mediator_ideal_equilibrium(mass,T)
    b=independent_bessel_eos(mass/T)
    for name in ("rho_over_T4","P_over_T4"):
        assert math.isclose(row[name],b[name],rel_tol=2e-8)
    assert 0<row["delta_gstar_rho"]<3

def test_light_vector_thermal_energy_conditional_not_actual_occupation():
    a=mediator_ideal_equilibrium(.03,1)
    assert a["delta_gstar_rho"]>2.99
    assert mediator_ideal_equilibrium(3,1)["delta_gstar_rho"]<2

def test_two_vector_thresholds_differ_for_mother_and_daughter():
    low=gauge_thresholds(mV_over_me=.03)
    assert low["chi1_anti_chi1_to_2V_open_at_rest"]
    assert not low["chi2_anti_chi2_to_2V_open_at_rest"]
    assert not low["V_to_e_minus_e_plus_open_at_rest"]
    high=gauge_thresholds(mV_over_me=.3)
    assert not high["chi1_anti_chi1_to_2V_open_at_rest"]
    assert high["V_to_two_active_nu_flavors_open"]
    assert low["yukawa_chi1_to_chi2_phi_open"]
    assert math.isclose(low["contact_G_times_mV_squared"],9e-7)


def test_dark_vector_decay_channels_open_in_distinct_mass_windows():
    low=gauge_thresholds(mV_over_me=.03)
    mid=gauge_thresholds(mV_over_me=.1)
    high=gauge_thresholds(mV_over_me=.3)
    assert not low["V_to_chi1_pair_open_at_rest"]
    assert not low["V_to_chi2_pair_open_at_rest"]
    assert mid["V_to_chi2_pair_open_at_rest"]
    assert not mid["V_to_chi1_pair_open_at_rest"]
    assert high["V_to_chi1_pair_open_at_rest"] and high["V_to_chi2_pair_open_at_rest"]
    assert low["chi1_pair_width_over_gchi2qchi2_mV"]==0
    assert mid["chi2_pair_width_over_gchi2qchi2_mV"]>0
    assert math.isclose(dirac_vector_width_coefficient(10,0),1/(12*math.pi))
