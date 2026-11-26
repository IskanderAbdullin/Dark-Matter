"""Pinned independent BEST massive phi source protocol, no false BEC claims."""
from __future__ import annotations
import math
import sys
from pathlib import Path

import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from best_massive_daughter_reference import massive_protocol,execute_best,BEST_PIN


def test_fixed_rest_width_normalizes_matrix_element_against_exact_phase_space():
    p=massive_protocol()
    m=p["mass_chi1"]
    k=p["daughter_rest_momentum"]
    matrix=p["spin_summed_fixed_width_matrix_element_squared"]
    assert math.isclose(matrix/(p["g_chi1"]*8*math.pi*m*m/k),1.,rel_tol=1e-14)
    assert math.isclose(
        p["independent_analytic_phi_number_source"],
        .0164959084596,rel_tol=1e-10)
    assert math.isclose(
        p["independent_analytic_parent_energy_loss"],
        .03646818570688685,rel_tol=1e-9)
    assert p["verified_best_git_commit_expected"]==BEST_PIN
    assert p["independent_exact_q0_phi_occupation_rate"]>10000
    assert p["not_physical_86eV_BEC_or_CLASS"]


def test_massless_limit_matrix_amplitude_agrees_with_prior_best_benchmark():
    p=massive_protocol(mphi=0.)
    assert math.isclose(
        p["spin_summed_fixed_width_matrix_element_squared"],
        16*math.pi*.1*4,rel_tol=1e-14)
    assert p["independent_exact_q0_phi_occupation_rate"]==0.


@pytest.mark.parametrize("r",[.01,.1,.2])
def test_fixed_gamma_number_is_mass_independent_but_phi_energy_changes(r):
    p=massive_protocol(mphi=.1*r)
    old=massive_protocol(mphi=0.)
    assert math.isclose(
        p["independent_analytic_phi_number_source"],
        old["independent_analytic_phi_number_source"],rel_tol=1e-9)
    assert p["spin_summed_fixed_width_matrix_element_squared"]>\
           old["spin_summed_fixed_width_matrix_element_squared"]
    assert p["independent_analytic_phi_energy_rate"]>\
           old["independent_analytic_phi_energy_rate"]


def test_bad_mass_threshold_and_degeneracies_rejected():
    with pytest.raises(ValueError):
        massive_protocol(mphi=.2)
    with pytest.raises(ValueError):
        massive_protocol(qmin=0)
    with pytest.raises(ValueError):
        massive_protocol(vegas_neval=0)


def test_missing_external_checkout_does_not_fake_result(tmp_path):
    with pytest.raises(ValueError,match="missing"):
        execute_best(massive_protocol(),tmp_path,tmp_path/"fake.json")
    assert not (tmp_path/"fake.json").exists()
