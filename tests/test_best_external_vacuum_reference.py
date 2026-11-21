"""BEST reference protocol is pinned and never mistakes absent run for evidence."""
from __future__ import annotations
import math
import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from best_external_vacuum_reference import benchmark_protocol,BEST_PIN,execute_best


def test_analytic_on_shell_vacuum_target_and_external_provenance():
    p=benchmark_protocol()
    assert p["external_run_verified"] is False
    assert "PROTOCOL ONLY" in p["status"]
    assert p["independent_solver"].startswith("BEST")
    assert math.isclose(
        p["independently_derived_vacuum_source_n_over_T3_per_time"],
        .0164959084596,rel_tol=1e-10)
    assert math.isclose(p["best_spin_summed_matrix_element_squared"],
                        16*math.pi*.1*1.*4,rel_tol=1e-15)
    assert BEST_PIN in p["best_source_url"]
    assert len(p["validation_requirements"])>=6


@pytest.mark.parametrize("bad",[
    {"mass_over_T":0.},{"parent_dof":0},{"rest_width":0.},
    {"grid_n":1},{"qmin":0.},{"qmax":1e-5},{"vegas_neval":10}])
def test_bad_model_inputs_rejected(bad):
    with pytest.raises(ValueError):
        benchmark_protocol(**bad)


def test_instantaneous_rate_scales_linearly_with_width():
    a=benchmark_protocol(rest_width=1.)
    b=benchmark_protocol(rest_width=2.)
    assert math.isclose(
        b["independently_derived_vacuum_source_n_over_T3_per_time"],
        2*a["independently_derived_vacuum_source_n_over_T3_per_time"],
        rel_tol=1e-15)


def test_missing_best_dependency_cannot_manufacture_run(tmp_path):
    with pytest.raises(ValueError,match="not found"):
        execute_best(benchmark_protocol(),tmp_path,tmp_path/"fake.json")
    assert not (tmp_path/"fake.json").exists()


def test_unversioned_best_source_is_rejected_before_import(tmp_path):
    (tmp_path/"besthep.py").write_text("# mock only; not a Git checkout")
    with pytest.raises(ValueError,match="provenance unavailable"):
        execute_best(benchmark_protocol(),tmp_path,tmp_path/"fake.json")
    assert not (tmp_path/"fake.json").exists()


def test_independent_parent_energy_flux_is_mass_times_width_times_number():
    from frw_reservoir_budget import thermal_mother_number
    p=benchmark_protocol()
    e=p["independently_derived_total_daughter_energy_gain_over_T4_per_time"]
    assert math.isclose(e,.1*thermal_mother_number(.1,g1=4),rel_tol=1e-13)
    assert math.isclose(e,.03646818570688685,rel_tol=1e-9)
    assert math.isclose(
        p["expected_parent_mean_energy_per_decay"],
        e/p["independently_derived_vacuum_source_n_over_T3_per_time"],
        rel_tol=1e-15)


def test_analytic_energy_loss_scales_with_width_and_mother_degeneracy():
    a=benchmark_protocol()
    b=benchmark_protocol(rest_width=2.)
    assert math.isclose(
        b["independently_derived_total_daughter_energy_gain_over_T4_per_time"],
        2*a["independently_derived_total_daughter_energy_gain_over_T4_per_time"],
        rel_tol=1e-15)
