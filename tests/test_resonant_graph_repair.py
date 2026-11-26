"""Independent rank/thermal-moment gate for repaired adaptive resonant meshes."""
from __future__ import annotations
import sys
from pathlib import Path

import pytest
from scipy.special import expit

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from resonant_graph_repair import (
    bridged_indices, exact_mother_thermal_moments,
    parent_grid_moments, evaluate_mesh
)
from collision_invariant_rank import audit_integer_grid,three_block_indices


@pytest.mark.parametrize("maximum",[5000,7000,8500,10000])
def test_bridging_cures_all_extra_graph_invariants(maximum):
    r=evaluate_mesh(max_index=maximum)
    assert r["graph_is_clean"],r
    assert r["collision_invariant_nullity"]==3
    assert r["extra_collision_invariants"]==0
    assert r["isolated_daughters"]==0
    assert r["discrete_mother_moments"]["parent_shells_all_positive"]


def test_adaptive_mesh_repair_is_not_equivalent_to_narrow_three_blocks():
    defective=audit_integer_grid(three_block_indices())
    fixed=audit_integer_grid(bridged_indices(max_index=5000))
    assert defective["active_graph_collision_invariants"]==4
    assert fixed["active_graph_collision_invariants"]==3
    assert fixed["active_unique_mother_energy_modes"]>5000
    assert fixed["daughter_grid_n"]<250


def test_static_parent_thermal_moments_converge_independently():
    ref=exact_mother_thermal_moments()
    results=[evaluate_mesh(max_index=m,reference=ref) for m in
             (5000,7000,8500,10000)]
    for key in ("number_over_T3","energy_over_T4","pressure_over_T4"):
        err=[abs(x["mother_relative_errors"][key]) for x in results]
        assert all(a>b for a,b in zip(err,err[1:])),(key,err)
    assert results[0]["parent_moment_max_abs_relative_error"]>.01
    assert results[2]["parent_moment_max_abs_relative_error"]<5e-5
    assert results[3]["parent_moment_max_abs_relative_error"]<3e-6
    assert abs(ref["number_over_T3"]-.3646818570688685)<1e-11


def test_no_condensate_or_class_claim_follows_from_good_graph():
    r=evaluate_mesh(max_index=8500)
    assert r["numerical_graph_and_moments_ready_for_further_studies"]
    assert r["NOT_a_validated_cosmological_PSD"]
    assert "No demonstrated full nonlinear evolution" in r["remaining"]


def test_input_validation_and_ir_support():
    with pytest.raises(ValueError):
        bridged_indices(core=1)
    with pytest.raises(TypeError):
        bridged_indices(core=2.5)
    x=bridged_indices(max_index=5000)
    mom=parent_grid_moments(x)
    assert mom["minimum_daughter_energy"]==.0005
    assert mom["maximum_daughter_energy"]==5.0005
    assert mom["parent_energy_min"]>=.1
