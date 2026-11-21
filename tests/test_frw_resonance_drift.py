"""Exact resonance failure when static 3-species energy bins are reused in FRW."""
from __future__ import annotations
import sys
from pathlib import Path
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from frw_resonance_drift import resonance_drift,run


def test_reference_resonance_is_exact_at_single_mass_temperature_ratio():
    r=resonance_drift(x0=.1,x=.1)
    assert r["number_of_physical_collision_invariants_at_x0"]==3
    assert r["maximum_absolute_energy_defect_per_reference_channel"]<1e-12
    assert r["rate_weighted_fraction_newly_forbidden_by_instantaneous_kinematics"]==0
    assert r["number_of_original_reactions_still_exactly_energy_resonant"]==r[
        "closed_graph_active_collisions"]


@pytest.mark.parametrize("new_x",[.11,.2,.5,1.,2.])
def test_no_fixed_massive_parent_comoving_momentum_remains_resonant(new_x):
    r=resonance_drift(x0=.1,x=new_x)
    assert r["minimum_signed_energy_defect"]>0
    assert r["maximum_signed_energy_defect"]>r["minimum_signed_energy_defect"]
    assert r["signed_weighted_energy_defect_relative_to_reference_source"]>0
    assert r["number_of_original_reactions_still_exactly_energy_resonant"]==0
    assert r["maximum_algebraic_energy_identity_residual"]<1e-12
    assert r["not_a_correct_FRW_collision_operator"]


def test_earlier_x_gives_negative_energy_defect_for_same_mother_momentum():
    r=resonance_drift(x0=.1,x=.05)
    assert r["maximum_signed_energy_defect"]<0
    assert r["signed_weighted_energy_defect_relative_to_reference_source"]<0
    assert r["number_of_original_reactions_still_exactly_energy_resonant"]==0


def test_severe_redshift_also_changes_allowed_reaction_graph():
    earlier=resonance_drift(x0=.1,x=.2)
    later=resonance_drift(x0=.1,x=2.)
    assert 0<=earlier["rate_weighted_fraction_newly_forbidden_by_instantaneous_kinematics"]<=1
    assert 0<later["rate_weighted_fraction_newly_forbidden_by_instantaneous_kinematics"]<=1
    assert later["rate_weighted_fraction_newly_forbidden_by_instantaneous_kinematics"]>=earlier[
        "rate_weighted_fraction_newly_forbidden_by_instantaneous_kinematics"]
    assert later["signed_weighted_energy_defect_relative_to_reference_source"]>earlier[
        "signed_weighted_energy_defect_relative_to_reference_source"]


def test_result_cannot_be_confused_with_actual_FRW_solution():
    report=run()
    assert len(report["cases"])==6
    assert "not a primordial spectrum" in report["caveat"]
    assert "Time-dependent on-shell remapping is mandatory" in report["physical_conclusion"]


@pytest.mark.parametrize("x0,x",[(0,.2),(.1,0),(-.1,.2)])
def test_bad_temperature_ratios_rejected(x0,x):
    with pytest.raises(ValueError):
        resonance_drift(x0=x0,x=x)
