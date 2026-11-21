"""Physical on-shell small-q vacuum source suppression and applicability."""
from __future__ import annotations
import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from vacuum_ir_tail_bound import thermal_vacuum_ir_fraction,run


@pytest.mark.parametrize("epsilon",[1e-5,1e-4,3e-4,5e-4,.001,.002,.01])
def test_rigorous_ir_vacuum_bound_holds(epsilon):
    result=thermal_vacuum_ir_fraction(epsilon)
    actual=result["actual_vacuum_fraction_below_threshold"]
    sharp=result["sharper_exponential_integral_bound"]
    simple=result["elementary_exponential_bound"]
    assert 0<=actual<=sharp*(1+1e-10)+1e-200
    assert sharp<=simple*(1+1e-10)+1e-200
    assert not result["applicable_to_stimulated_nonlinear_source"]


def test_numeric_low_q_fraction_from_independent_1d_quad():
    result=thermal_vacuum_ir_fraction(.0005)
    assert math.isclose(result["actual_vacuum_fraction_below_threshold"],
                        6.106102547396695e-7,rel_tol=2e-5)
    result2=thermal_vacuum_ir_fraction(.001)
    assert math.isclose(result2["actual_vacuum_fraction_below_threshold"],
                        2.37124370929558e-5,rel_tol=2e-5)


def test_total_source_matches_independent_on_shell_vacuum_rate():
    result=thermal_vacuum_ir_fraction(.0005)
    assert math.isclose(result["total_vacuum_source_per_T3_Gamma"],
                        .0164959084596,rel_tol=1e-10)


def test_monotonicity_and_empty_daughter_only_warning():
    result=run()
    fractions=[x["actual_vacuum_fraction_below_threshold"] for x in result["cases"]]
    assert fractions==sorted(fractions)
    assert all(x["applicable_to_stimulated_nonlinear_source"] is False
               for x in result["cases"])
    assert result["not_a_cosmological_freezeout"]


def test_bad_source_parameters_rejected():
    for epsilon,m in [(0.,.1),(-.1,.1),(.001,-.1)]:
        with pytest.raises(ValueError):
            thermal_vacuum_ir_fraction(epsilon,mass=m)
