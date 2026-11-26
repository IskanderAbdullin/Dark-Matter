"""Conditional SM/dark separate comoving entropy, not a kinetic prediction."""
import math
import pytest
from kinetics.phi_bose_dark_sm_entropy_split_gate import (
    relativistic_dark_gstar_s,decoupled_entropy_temperature_ratio,benchmark,
)


def test_initial_dark_five_population_relativistic_entropy_gstar_exact():
    assert relativistic_dark_gstar_s()==8.
    assert decoupled_entropy_temperature_ratio(10.75,10.75,8,8)==1.


def test_independent_comoving_entropy_conservation_between_sectors():
    # With a shared scale factor growth a and no comoving entropy exchange,
    # gsm(Tsm)*Tsm³*a³ and gd(Td)*Td³*a³ each remain constant.
    a_growth=50.
    g0SM=61.75;g1SM=10.75
    g0D=8.;g1D=1.
    T0=1.
    Tsm=T0/a_growth*(g0SM/g1SM)**(1/3)
    Td=T0/a_growth*(g0D/g1D)**(1/3)
    assert math.isclose(Td/Tsm,
                        decoupled_entropy_temperature_ratio(
                            g0SM,g1SM,g0D,g1D),rel_tol=2e-14)
    assert math.isclose(g0SM*T0**3,
                        g1SM*Tsm**3*a_growth**3,rel_tol=2e-14)
    assert math.isclose(g0D*T0**3,
                        g1D*Td**3*a_growth**3,rel_tol=2e-14)


def test_conditional_extreme_temperature_ratios_are_not_an_actual_Dirac_relic_history():
    b=benchmark()
    assert math.isclose(b["no_dark_internal_entropy_release_xi"],
                        .5583722596874463,rel_tol=2e-14)
    assert math.isclose(b["complete_hypothetical_dark_reheating_into_phi_xi"],
                        1.1167445193748926,rel_tol=2e-14)
    assert b["ratio_of_two_cases"]==pytest.approx(2.,rel=3e-15)
    assert not b["actual_chi2_relic_compatible_with_gD_after_1"]
    assert not b["class_ready"]


def test_reject_malformed_entropy_DoF():
    with pytest.raises(ValueError):
        decoupled_entropy_temperature_ratio(gd_after=0)
    with pytest.raises(ValueError):
        decoupled_entropy_temperature_ratio(gsm_before=float("nan"))


def test_finite_mass_FD_entropy_integral_quantifies_true_dark_gstar_near_crossing():
    from kinetics.phi_bose_dark_sm_entropy_split_gate import (
        physical_FD_dirac_entropy_gstar_s,physical_initial_dark_gstar_s,
    )
    assert math.isclose(physical_FD_dirac_entropy_gstar_s(0,4.),
                        3.5,rel_tol=2e-12)
    assert physical_FD_dirac_entropy_gstar_s(.15,4.)<3.5
    assert physical_FD_dirac_entropy_gstar_s(10,4.)<1.
    t=physical_initial_dark_gstar_s()
    assert math.isclose(t,7.991670794968991,rel_tol=2e-10)
    assert 0<8./t-1<.002
    b=benchmark()
    assert math.isclose(
       b["actual_thermal_dark_gstar_s_with_mass_at_transport_crossing_T_0p688GeV"],t,
       rel_tol=1e-12)
    assert b["initial_massless_approx_error_fraction"]<.002
