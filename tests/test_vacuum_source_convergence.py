"""Independent continuum rate, finite-cell ODE initial flux and nonlinear error."""
from __future__ import annotations
import math
import sys
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from vacuum_source_convergence import one_source,run


def test_initial_closed_rhs_equals_sum_of_on_shell_thermal_pair_rates():
    for max_index in (5000,8500,10000):
        a=one_source(max_index=max_index)
        assert a["reaction_graph_invariants"]==3
        assert abs(a["ODE_initial_phi_number_rate"]-
                   a["initial_grid_source_rate"])<1e-11
        assert math.isclose(a["continuum_analytic_source_rate"],
                            .0164959084596,rel_tol=1e-10)


def test_uv_source_error_decreases_over_selected_tested_meshes():
    x=[one_source(max_index=m) for m in (5000,7000,8500,10000)]
    err=[abs(z["grid_source_relative_error"]) for z in x]
    assert all(a>b for a,b in zip(err,err[1:])),err
    assert .015<err[0]<.018
    assert err[-1]<.0003


def test_static_mother_number_accuracy_is_not_source_rate_accuracy():
    r=one_source(max_index=8500)
    fd=abs(r["initial_parent_number_relative_error"])
    src=abs(r["grid_source_relative_error"])
    assert fd<1e-5
    assert src>fd*50
    assert src>.0001


def test_more_daughter_bins_do_not_necessarily_improve_source_error():
    a=one_source(core=100,stride=50,max_index=10000)
    b=one_source(core=100,stride=25,max_index=10000)
    assert b["daughter_modes"]>a["daughter_modes"]
    assert abs(b["grid_source_relative_error"])>abs(a["grid_source_relative_error"])


def test_early_nonlinear_production_can_be_more_cutoff_sensitive_than_vacuum():
    a=run(quick=True,dynamic=True)
    v=a["optional_cutoff_sensitivity_comparison"]
    assert v is not None
    assert v["initial_vacuum_rate_fractional_change"]>0
    assert v["nonlinear_tGamma0p001_phi_number_fractional_change"]>0
    assert v["ratio_of_two_cutoff_sensitivities"]>3
    assert "Not a physical Bose enhancement factor" in v["interpretation"]
