"""Conservation, cutoff and BDF-tolerance discrimination for closed kinetics."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from cutoff_and_time_convergence import run,one_run,relative_delta


def test_reproducible_short_closed_source_with_physical_mother_mass():
    result=one_run(5000)
    assert result["graph_invariants"]==3
    assert result["daughter_modes"]==199
    assert result["mother_modes"]==5099
    assert .008<result["dynamic_phi_number_over_T3"]<.0095
    assert .003<result["dynamic_phi_mean_q_over_T"]<.005
    assert result["maximum_absolute_charge_or_energy_drift"]<1e-8
    assert abs(result["initial_parent_FD_relative_errors"]["number_over_T3"])<.004


def test_nonconverged_upper_momentum_cutoff_dominates_time_discretization():
    result=run(quick=True)
    a,b=result["cutoff_variations"]
    assert a["max_index"]==5000
    assert b["max_index"]==7000
    assert b["dynamic_phi_number_over_T3"]>a["dynamic_phi_number_over_T3"]
    assert result["relative_number_cutoff_discrepancy_last_two"]>.05
    assert result["relative_number_time_discrepancy"]<.003
    assert result["ratio_cutoff_discrepancy_to_time_discrepancy"]>10
    assert result["initial_FD_moments_are_not_proxy_for_dynamic_source_accuracy"]


def test_cutoff_relative_delta_sign_and_normalization():
    assert relative_delta(11,10)==.1 or abs(relative_delta(11,10)-.1)<1e-12
    assert relative_delta(3,3)==0
