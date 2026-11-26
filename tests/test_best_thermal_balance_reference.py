"""Exact Fermi/Bose thermal detailed balance, independent of finite-grid ODE."""
from __future__ import annotations
import math
import sys
from pathlib import Path
import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from best_thermal_balance_reference import fermi,bose,event_balance,protocol,execute_best


@pytest.mark.parametrize("qphi,qchi",[
    (.005,.6),(.01,.3),(.03,.2),(.1,.1),(.5,.1),(1.,1.),(3.,2.)])
def test_exact_microscopic_quantum_forward_reverse_thermal_balance(qphi,qchi):
    r=event_balance(qphi,qchi)
    assert r["chemical_equilibrium"] is True
    assert abs(r["relative_balance_error"])<1e-12
    assert r["forward_quantum_factor"]>0
    assert r["inverse_quantum_factor"]>0


def test_wrong_parent_chemical_potential_gives_nonzero_signed_collision():
    a=event_balance(.1,.2)
    b=event_balance(.1,.2,mu_parent=-.1)
    c=event_balance(.1,.2,mu_parent=-.3)
    assert abs(a["relative_balance_error"])<1e-12
    assert b["relative_balance_error"]>.01
    assert c["relative_balance_error"]<-.01


def test_analytic_thermal_protocol_is_unrun_and_uses_finite_bose_IR():
    p=protocol()
    assert "NOT EXECUTED" in p["status"]
    assert not p["external_best_executed"]
    assert p["not_cosmological_freezeout"]
    assert p["mu_parent"]==p["mu_phi"]+p["mu_chi2"]
    assert p["mu_phi"]<0
    assert p["max_pointwise_analytic_gain_loss_error"]<1e-11
    assert p["number_of_independent_analytic_events"]==40
    assert p["gphi"]==1 and p["gchi2"]==4 and p["g1"]==4


def test_fermion_bounds_and_positive_bose_occupation_for_safe_mu():
    energies=np.geomspace(.001,10,100)
    f=fermi(energies,.1)
    b=bose(energies,-.3)
    assert np.all((f>0)&(f<1))
    assert np.all(b>0)
    assert np.all(np.diff(b)<0)


def test_invalid_bose_chemical_potential_or_offshell_event_rejected():
    with pytest.raises(ValueError):
        bose(.1,mu=.2)
    with pytest.raises(ValueError):
        event_balance(.001,.001)
    with pytest.raises(ValueError):
        event_balance(-.1,1.)


def test_missing_external_best_checkout_is_not_silently_replaced(tmp_path):
    with pytest.raises(ValueError,match="not found"):
        execute_best(tmp_path,tmp_path/"output.json")
    assert not (tmp_path/"output.json").exists()
