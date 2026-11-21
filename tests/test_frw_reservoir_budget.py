"""Necessary parent-number budget under an explicitly CLOSED interpretation."""
from __future__ import annotations

import math
import sys
from pathlib import Path
from scipy.special import zeta

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from frw_reservoir_budget import mother_budget, thermal_mother_number
from frw_fixed_bath_pilot import FRWFixedBathPilot


def test_thermal_parent_abundance_massless_exact_limit():
    answer=thermal_mother_number(0.)
    assert math.isclose(answer,3*zeta(3,1)/math.pi**2,rel_tol=1e-12)


def test_thermal_parent_density_decreases_with_mass_over_temperature():
    vals=[thermal_mother_number(x) for x in (0.,.05,.1,1.,5.)]
    assert all(a>b>0 for a,b in zip(vals,vals[1:]))
    assert math.isclose(vals[1],.3652061032773952,rel_tol=1e-10)


def test_closed_daughter_charge_implies_parent_replenishment_lower_bound():
    parent=thermal_mother_number(.05)
    low=mother_budget(.1)
    impossible=mother_budget(parent+.01)
    assert low["closed_without_external_parent_source_possible_by_number"]
    assert math.isclose(low["closed_reaction_minimum_parent_depletion_fraction"],.1/parent,rel_tol=1e-12)
    assert not impossible["closed_without_external_parent_source_possible_by_number"]
    assert math.isclose(impossible["minimum_required_external_net_parent_replenishment_over_T3"],.01,abs_tol=1e-12)


def test_strong_fixed_bath_frw_exceeds_closed_initial_mother_number():
    # Valid fixed-bath calculation, invalid to relabel as CLOSED.
    result=FRWFixedBathPilot(nq=80,kappa=.2).integrate()
    diag=mother_budget(result["final"]["nphi_over_T3"],.05)
    assert not diag["closed_without_external_parent_source_possible_by_number"]
    assert diag["minimum_required_external_net_parent_replenishment_over_T3"]>.002
    assert result["frozen_number_to_1e-4_over_late_window"]
