"""Mother resonant mode-set changes handled by conservative q-remapping, or fail closed."""
from __future__ import annotations
import sys
from pathlib import Path
import math
import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from narrow_closed_frw_split import NarrowClosedFRWPilot


def test_two_parent_modes_disappear_but_no_physical_mother_number_disappears():
    m=NarrowClosedFRWPilot(x0=.02,xend=.022,steps=4,kappa=10.,
                           allow_mother_mode_changes=True)
    assert m.models[0].nparent==385
    assert m.models[-1].nparent==383
    assert all(mod.graph["active_graph_collision_invariants"]==3
               for mod in m.models)
    assert all(m.parent_edges[j][-1]==m.qmax
               for j in range(len(m.parent_edges)))
    outcome=m.run()
    assert outcome["parent_mode_changes_permitted"]
    assert outcome["final_shared_parent_modes"]<outcome["shared_parent_modes"]
    assert any(row["parent_mode_count_change"]<0 for row in outcome["steps"])
    assert outcome["final"]["n_phi"]>0
    assert outcome["max_global_charge_drift"]<1e-8
    for row in outcome["steps"]:
        assert row["maximum_charge_drift"]<1e-8
        assert abs(row["collision_discrete_energy_error"])<1e-8
        assert math.isfinite(row["finite_volume_remap_energy_error"])
        assert math.isfinite(row["collisionless_FRW_energy_change"])
        assert row["min_fermion_occupation"]>=-1e-8
        assert row["max_fermion_occupation"]<=1+1e-8


def test_mode_changes_are_rejected_without_explicit_mode_switch():
    with pytest.raises(ValueError,match="energies changed"):
        NarrowClosedFRWPilot(x0=.02,xend=.022,steps=4,kappa=10.,
                              allow_mother_mode_changes=False)


def test_rank_gate_blocks_larger_x_when_graph_has_artificial_invariants():
    with pytest.raises(ValueError,match="reaction invariants"):
        NarrowClosedFRWPilot(x0=.02,xend=.03,steps=4,kappa=10.,
                              allow_mother_mode_changes=True)


def test_each_instantaneous_parent_energy_is_physical_and_on_shell():
    model=NarrowClosedFRWPilot(x0=.02,xend=.022,steps=4,kappa=10.,
                               allow_mother_mode_changes=True)
    for m,x in zip(model.models,model.xgrid):
        assert m.parent_E[0]>x
        assert m.max_resonance_error<1e-12
        assert m.graph["active_graph_collision_invariants"]==3
