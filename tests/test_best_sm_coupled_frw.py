"""Pinned BEST + finite-SM-heat coupled FRW: fail-closed execution protocol."""
from __future__ import annotations
import math
import sys
from pathlib import Path
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from best_sm_coupled_frw import sm_coupled_protocol,execute_best
from best_external_vacuum_reference import BEST_PIN


def test_true_dynamic_three_species_plus_sm_heat_budget_specification():
    p=sm_coupled_protocol()
    assert p["kinetic_source"]["verified_best_git_commit_expected"]==BEST_PIN
    assert p["kinetic_source"]["mass_chi1"]==.1
    assert p["kinetic_source"]["mass_phi"]==.001
    assert p["integrator"]=="exprb2"
    assert p["steps"]==8
    assert math.isclose(p["end_tGamma"],4e-5,rel_tol=1e-14)
    assert math.isclose(p["rd_final_expected_a"],math.sqrt(1.4),rel_tol=1e-14)
    assert p["elastic_gamma_el_over_Gamma_rest"]==1e5
    assert p["SM_radiation_gstarS_equals_gstarRho"]==55.
    assert p["SM_energy_and_temperature_dynamic"]
    assert not p["real_dark_matter_86eV_or_CLASS_relic_validated"]
    assert not p["external_BEST_executed"]


def test_elastic_contact_is_number_conserving_not_hidden_mother_creation():
    p=sm_coupled_protocol(gamma_el=0.)
    assert p["elastic_gamma_el_over_Gamma_rest"]==0
    assert p["SM_thermalization_is_number_conserving_BGK_not_microphysical_2to2"]
    assert p["FRW_Hubble_prescribed_not_sourced_by_SM_or_dark"]
    assert "elastic half-step" in p["time_splitting"]


@pytest.mark.parametrize("params",[
    {"gamma_el":-1},{"gamma_el":float("nan")},
    {"gstar":0},{"gstar":float("inf")},
    {"t0":0},{"steps":0},{"dtGamma":0},
    {"nq":1},{"neval":100},{"method":"unsupported"},
])
def test_invalid_or_falsely_physical_SM_BEST_protocol_is_rejected(params):
    with pytest.raises(ValueError):
        sm_coupled_protocol(**params)


def test_external_BEST_cannot_be_spoofed_by_a_missing_directory(tmp_path):
    output=tmp_path/"fake_coupled_SM.json"
    with pytest.raises(ValueError,match="missing"):
        execute_best(tmp_path,output,sm_coupled_protocol(steps=1))
    assert not output.exists()
