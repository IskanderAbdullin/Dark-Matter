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
