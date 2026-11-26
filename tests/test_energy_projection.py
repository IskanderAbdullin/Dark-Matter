"""Detect energy projection inconsistency in occupation-dependent subevents."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from energy_projection_audit import diagnose


def test_vacuum_events_project_energy_exactly():
    a=diagnose(50,64,48,populated=False)
    assert abs(a["signed_relative_energy_discrepancy"]) < 1e-11
    assert a["absolute_per_parent_event_relative_error"] < 1e-11


def test_occupation_dependent_rates_break_event_energy_conservation():
    coarse=diagnose(50,64,48,populated=True)
    fine=diagnose(100,140,72,populated=True)
    assert abs(coarse["signed_relative_energy_discrepancy"]) > 0.005
    assert 0 < fine["absolute_per_parent_event_relative_error"] < coarse["absolute_per_parent_event_relative_error"]
