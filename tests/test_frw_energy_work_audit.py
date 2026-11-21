"""FRW pressure work versus remap and collision physical-energy defects."""
from __future__ import annotations
import sys
from pathlib import Path
import math
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from frw_energy_work_audit import (
    physical_parent_pressure,direct_frw_pressure_work,verify_free_expansion,
    audit_split_ledger,run
)
from parent_q_remap import thermal_shell_averages
from narrow_closed_frw_split import NarrowClosedFRWPilot
from strang_frw_split import SymmetricFRWPilot


@pytest.mark.parametrize("x",[.05,.1,1.,5.,20.])
def test_exact_FV_massive_parent_FRW_pressure_identity(x):
    q=np.r_[0.,np.geomspace(.01,20.,80)]
    f=thermal_shell_averages(q,.05)
    result=physical_parent_pressure(q,f,x)
    assert abs(result["pressure_continuity_residual"])<1e-11
    assert result["rho_over_T4"]>0
    assert 0<result["pressure_over_T4"]<result["rho_over_T4"]/3


@pytest.mark.parametrize("x1",[.1,.5,2.,10.])
def test_independent_integral_over_expansion_equals_energy_change(x1):
    q=np.r_[0.,np.geomspace(.01,30.,120)]
    f=thermal_shell_averages(q,.05)
    report=verify_free_expansion(q,f,.05,x1)
    assert abs(report["independent_integration_residual"])<3e-10
    assert report["integral_of_rho_minus_3P_over_x"]>0


def test_lie_split_physical_energy_ledger_is_complete_not_energy_closed():
    sim=NarrowClosedFRWPilot(steps=2).run()
    result=audit_split_ledger(sim)
    assert abs(result["ledger_residual"])<1e-10
    assert result["maximum_individual_ledger_residual"]<1e-10
    assert result["unphysical_energy_error_is_not_assumed_zero"]
    assert result["not_a_self_consistent_cosmology"]
    assert len(sim["steps"])==2


def test_symmetric_split_physical_energy_ledger_is_complete_not_energy_closed():
    sim=SymmetricFRWPilot(steps=2).run()
    result=audit_split_ledger(sim)
    assert abs(result["ledger_residual"])<1e-10
    assert result["maximum_individual_ledger_residual"]<1e-10
    assert result["unphysical_energy_error_is_not_assumed_zero"]
    assert result["not_a_self_consistent_cosmology"]


def test_pressure_work_zero_for_unchanged_x():
    q=[0.,.1,1.,2.]
    assert direct_frw_pressure_work(q,[.1,.2,.3],.2,.2)==0.


def test_invalid_inputs_and_empty_ledgers_are_rejected():
    with pytest.raises(ValueError):
        audit_split_ledger({"steps":[]})
    with pytest.raises(ValueError):
        direct_frw_pressure_work([0,1],[.5],2.,1.)
    with pytest.raises(ValueError):
        physical_parent_pressure([0,1],[.5],-1.)


def test_cli_report_never_claims_full_frw_energy_conservation():
    report=run()
    assert len(report["cases"])==4
    assert "not proof" in report["warning"]
    assert all(abs(x["independent_integration_residual"])<3e-10
               for x in report["cases"])
