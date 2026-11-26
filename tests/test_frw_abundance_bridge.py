"""Relic-amplitude crosscheck against the existing exact thermal Bose benchmark."""
from __future__ import annotations
import math
import sys
from pathlib import Path

import pytest
from scipy.special import zeta

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"cosmology"))
from frw_abundance_bridge import relic_from_freezeout_number, required_freezeout_number, N_BE
from thermal_bose_class import thermal_background


def test_thermal_limit_exactly_matches_existing_verified_cosmology_anchor():
    baseline=thermal_background(m=86,gdec=55,gtoday=3.909)
    got=relic_from_freezeout_number(zeta(3,1)/math.pi**2)
    assert math.isclose(got["hypothetical_present_Omega_phi_h2"],
                        baseline["omega_h2"],rel_tol=1e-14)
    assert math.isclose(N_BE,zeta(3,1)/math.pi**2,rel_tol=1e-15)


def test_mass_and_entropy_scalings():
    default=relic_from_freezeout_number(.05)
    double_mass=relic_from_freezeout_number(.05,m_eV=172)
    double_gstar=relic_from_freezeout_number(.05,gstarS_dec=110)
    assert math.isclose(double_mass["hypothetical_present_Omega_phi_h2"],
                        2*default["hypothetical_present_Omega_phi_h2"],rel_tol=1e-13)
    assert math.isclose(double_gstar["hypothetical_present_Omega_phi_h2"],
                        .5*default["hypothetical_present_Omega_phi_h2"],rel_tol=1e-13)


def test_inverse_required_number_gives_target():
    q=required_freezeout_number(target_omega=.12)
    d=relic_from_freezeout_number(q)
    assert .12<q<.13
    assert math.isclose(d["hypothetical_present_Omega_phi_h2"],.12,rel_tol=1e-14)
    assert 0.33<q/.3652061032773952<.35


def test_no_shape_or_neff_claim_from_relic_number_only():
    d=relic_from_freezeout_number(.10)
    assert "not_inferred" in d
    assert not any(key.lower().startswith("deltaneff") for key in d)


@pytest.mark.parametrize("n",[-1,.0])
def test_nonnegative_moments(n):
    if n<0:
        with pytest.raises(ValueError):relic_from_freezeout_number(n)
    else:
        assert relic_from_freezeout_number(n)["hypothetical_present_Omega_phi_h2"]==0
