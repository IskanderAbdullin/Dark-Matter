"""Finite-volume massive-mother FRW shell quadrature and continuity regression."""
from __future__ import annotations
import math
import sys
from pathlib import Path
import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from frw_mother_finite_volume import CollisionlessMotherFV,run
from frw_collisionless_mother import frozen_mother_background


@pytest.mark.parametrize("x",[.05,.1,1.,2.,5.,20.])
def test_fixed_q_bins_conserve_number_and_kinetic_energy_continuity(x):
    model=CollisionlessMotherFV(nq=100)
    result=model.moments(x)
    assert math.isclose(result["n_over_T3"],model.number,rel_tol=1e-14)
    assert abs(result["x_d_rho_dx_minus_rho_plus_3P"])<1e-10
    assert abs(result["d_ln_physical_rho_by_d_ln_a"]-
               result["expected_FRW_energy_slope"])<1e-10
    assert 0<result["w"]<1/3
    assert np.all(model.f_i<=1) and np.all(model.f_i>=0)


def test_number_fv_matches_exact_infinite_FD_integral_on_even_coarse_grids():
    ref=frozen_mother_background(5.,x_dec=.05)
    for n in (30,50,100,200,400):
        model=CollisionlessMotherFV(nq=n)
        assert abs(model.number/ref["comoving_mother_number_over_T3"]-1)<1e-10


def test_conserving_number_does_not_imply_converged_pressure():
    errors=[]
    for n in (30,50,100,200,400,600):
        m=CollisionlessMotherFV(nq=n)
        er=m.reference_error(5.)
        assert abs(er["n_relative_error"])<1e-10
        errors.append(er)
    rho=[abs(row["rho_relative_error"]) for row in errors]
    p=[abs(row["pressure_relative_error"]) for row in errors]
    assert all(a>b for a,b in zip(rho,rho[1:])),rho
    assert all(a>b for a,b in zip(p,p[1:])),p
    assert p[0]>.1
    assert p[-1]<.0005
    assert rho[-1]<.0001


def test_independent_finite_difference_matches_FV_drho_dx():
    model=CollisionlessMotherFV(nq=100)
    x=2.;h=1e-4
    derivative=(model.moments(x+h)["rho_over_T4"]-
                model.moments(x-h)["rho_over_T4"])/(2*h)
    assert math.isclose(derivative,model.moments(x)["analytic_d_rho_over_T4_dx"],
                        rel_tol=1e-9,abs_tol=1e-10)


def test_initial_occupations_are_not_rethermalized():
    model=CollisionlessMotherFV(nq=100)
    original=model.f_i.copy()
    for x in (.05,1.,5.,20.):
        _=model.moments(x)
    assert np.array_equal(original,model.f_i)


def test_cli_survey_separates_source_free_from_closed_kinetics():
    r=run()
    assert len(r["cases"])==6
    assert "SOURCE-FREE" in r["status"]
    assert "Combine" in r["next"]


def test_invalid_grid_or_x_parameters_rejected():
    with pytest.raises(ValueError):
        CollisionlessMotherFV(nq=1)
    with pytest.raises(ValueError):
        CollisionlessMotherFV(nq=50,qmin=-1)
    with pytest.raises(ValueError):
        CollisionlessMotherFV(nq=50).moments(0.)
