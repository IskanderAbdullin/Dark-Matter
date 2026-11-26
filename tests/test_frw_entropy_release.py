"""Entropy-aware collisionless FRW Liouville for massive decoupled parent."""
from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from frw_entropy_release import SmoothEntropyStep,frozen_entropy_mother,run
from frw_collisionless_mother import frozen_mother_background


@pytest.mark.parametrize("x",[.05,.2,.5,1.,2.,5.,10.,20.])
def test_nonconstant_SM_entropy_conserves_mother_comoving_number_and_energy(x):
    r=frozen_entropy_mother(x)
    n0=r["initial_mother_n_over_Tdec3"]
    assert math.isclose(r["comoving_number_charge"],n0,rel_tol=1e-11)
    assert math.isclose(r["current_mother_n_over_T3"]/n0,
                        r["gstarS_number_ratio"],rel_tol=1e-11)
    assert abs(r["entropy_aware_energy_identity_residual"])<1e-10
    assert abs(r["physical_d_ln_rho_d_ln_a"]-
               r["continuity_expected_slope"])<1e-10
    assert -4<r["continuity_expected_slope"]<-3


def test_constant_entropy_degrees_recover_previous_exact_liouville_benchmark():
    entropy=SmoothEntropyStep(g_early=55.,g_late=55.)
    for x in (.05,.5,5.):
        a=frozen_entropy_mother(x,entropy=entropy)
        b=frozen_mother_background(x)
        assert math.isclose(a["current_mother_n_over_T3"],
                            b["comoving_mother_number_over_T3"],rel_tol=1e-12)
        assert math.isclose(a["rho_over_T4"],
                            b["dimensionless_rho_over_T4"],rel_tol=1e-12)
        assert math.isclose(a["pressure_over_T4"],
                            b["dimensionless_pressure_over_T4"],rel_tol=1e-12)
        assert a["entropy_dln_gstar_dln_x"]==0
        assert a["p_over_T_to_comoving_q_scale"]==1.


def test_analytic_entropy_step_derivative_and_FRW_energy_via_finite_difference():
    entropy=SmoothEntropyStep()
    x=1.;h=1e-5
    gplus=entropy.value(x*math.exp(h))
    gminus=entropy.value(x*math.exp(-h))
    alphafinite=(math.log(gplus)-math.log(gminus))/(2*h)
    r=frozen_entropy_mother(x,entropy=entropy)
    assert math.isclose(alphafinite,r["entropy_dln_gstar_dln_x"],
                        rel_tol=1e-9)
    rhoplus=frozen_entropy_mother(x*math.exp(h))["rho_over_T4"]
    rhominus=frozen_entropy_mother(x*math.exp(-h))["rho_over_T4"]
    numerical_derivative=(rhoplus-rhominus)/(2*h)
    assert math.isclose(numerical_derivative,
                        r["x_drho_over_T4_dx_from_integral"],
                        rel_tol=1e-8)


def test_entropy_release_changes_n_over_T3_but_not_n_a3():
    at_dec=frozen_entropy_mother(.05)
    late=frozen_entropy_mother(20.)
    ratio=late["current_mother_n_over_T3"]/at_dec["current_mother_n_over_T3"]
    assert .070<ratio<.072
    assert late["p_over_T_to_comoving_q_scale"]<.42
    assert late["a_over_a_decoupling"]>900
    assert math.isclose(at_dec["comoving_number_charge"],
                        late["comoving_number_charge"],rel_tol=1e-11)


def test_invalid_degrees_or_temperature_rejected():
    with pytest.raises(ValueError):
        SmoothEntropyStep(g_early=2,g_late=4)
    with pytest.raises(ValueError):
        frozen_entropy_mother(0.)
    with pytest.raises(ValueError):
        SmoothEntropyStep().value(-1.)


def test_full_synthetic_entropy_history_is_not_class_model():
    r=run()
    assert len(r["cases"])==8
    assert r["not_a_physical_DM_relic_prediction"]
    assert "test function" in r["entropy_model"]
