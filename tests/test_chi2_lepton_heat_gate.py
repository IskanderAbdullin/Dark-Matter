"""Scoped χ2 vector heat gate, NOT thermal χ2 abundance prediction."""
import pytest
from kinetics.chi2_lepton_heat_gate import audit


def test_chi2_equal_charge_portal_has_comparable_thermal_heat():
    a=audit(nr=5,nmu=5)
    for row in a["moments"].values():
        assert row["Q1"]>0 and row["Q2"]>0
        assert .5 < row["Q2_over_Q1"] < 2.
    assert a["moments"]["two_active_nu"]["Q2_over_Q1"]>1.2


def test_chi2_decay_threshold_86ev_scalar():
    with pytest.raises(ValueError):
        audit(M1_me=.1,M2_me=.1,mphi_eV=86.)
