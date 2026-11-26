"""Changing M/T can create spurious collision invariants; repair by UV support."""
from __future__ import annotations
import sys
from pathlib import Path
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from dynamic_rank_uv_repair import repaired_instantaneous_grid,run


@pytest.mark.parametrize("mass,expected_upper",[
    (.02,300),(.025,400),(.03,500),(.04,1000),(.05,1500),(.08,3500),(.1,5000)
])
def test_minimum_known_uv_index_restores_three_invariants(mass,expected_upper):
    result=repaired_instantaneous_grid(mass)
    assert result["selected_upper_index"]==expected_upper
    assert result["number_of_collision_invariants"]==3
    assert result["collision_graph_pass"]
    assert result["attempts"][-1]["graph_pass"]
    assert not result["frozen_DM_spectrum_available"]


def test_smallest_clean_graph_has_bad_thermal_parent_normalization():
    result=repaired_instantaneous_grid(.025)
    assert not result["independent_parent_moment_accuracy_better_than_one_percent"]
    assert max(abs(v) for v in result["mother_thermal_relative_errors"].values())>.1


def test_rank_gate_rejects_insufficient_upper_cutoff_instead_of_guessing():
    with pytest.raises(ValueError,match="No invariant-clean"):
        repaired_instantaneous_grid(.1,upper_candidates=(300,1000,2000,3000))


def test_graph_repair_scan_lists_each_failed_attempt():
    result=repaired_instantaneous_grid(.03)
    assert [r["upper_index"] for r in result["attempts"]]==[300,400,500]
    assert result["attempts"][0]["extra_invariants"]>0
    assert result["attempts"][1]["extra_invariants"]>0
    assert result["attempts"][-1]["extra_invariants"]==0


def test_invalid_arguments_and_cosmology_gate():
    with pytest.raises(ValueError):
        repaired_instantaneous_grid(0)
    with pytest.raises(ValueError):
        repaired_instantaneous_grid(.1,upper_candidates=())
    report=run(quick=True)
    assert len(report["instantaneous_mass_sweep"])==4
    assert "no CLASS" in report["scope"]
