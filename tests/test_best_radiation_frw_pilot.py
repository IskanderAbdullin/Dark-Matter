"""Independent FRW physical pressure-work protocol is not a CLASS relic."""
from __future__ import annotations
import math
import sys
from pathlib import Path
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from best_radiation_frw_pilot import rd_protocol,execute_best
from best_external_vacuum_reference import BEST_PIN


def test_rd_scale_factor_and_mass_over_temperature_are_physically_consistent():
    p=rd_protocol()
    assert p["rd_t0"]==.0001
    assert p["rd_initial_scale_factor"]==1
    assert math.isclose(p["rd_final_expected_a"],math.sqrt(1.4),rel_tol=1e-14)
    assert p["rd_initial_H_over_Gamma"]==5000
    assert math.isclose(p["end_tGamma"],4e-5,rel_tol=1e-14)
    assert p["kinetic_source"]["verified_best_git_commit_expected"]==BEST_PIN
    assert p["kinetic_source"]["mass_chi1"]==.1
    assert p["kinetic_source"]["mass_phi"]==.001
    assert p["dynamic_mother_is_NOT_instantaneously_rethermalized"]
    assert not p["primordial_86eV_relic_or_CLASS_ready"]


def test_halving_dt_preserves_exact_prescribed_final_scale_factor():
    coarse=rd_protocol(steps=4,dtGamma=1e-5)
    fine=rd_protocol(steps=8,dtGamma=5e-6)
    assert math.isclose(coarse["rd_final_expected_a"],
                        fine["rd_final_expected_a"],rel_tol=1e-14)
    assert coarse["end_tGamma"]==fine["end_tGamma"]


def test_missing_best_checkout_cannot_forge_prescribed_frw_evolution(tmp_path):
    with pytest.raises(ValueError,match="missing"):
        execute_best(tmp_path,tmp_path/"fake.json",rd_protocol())
    assert not (tmp_path/"fake.json").exists()


@pytest.mark.parametrize("arg",[
    {"t0":0},{"steps":0},{"dtGamma":0},
    {"method":"bogus"},{"nq":1},{"neval":0},{"phi_mass":.2}
])
def test_invalid_rd_physical_inputs_fail_closed(arg):
    with pytest.raises(ValueError):
        rd_protocol(**arg)


def test_external_rd_is_conditional_not_conservation_proof():
    p=rd_protocol()
    assert "NOT EXECUTED" in p["status"]
    assert p["full_SM_energy_transfer_or_entropy_included"] is False
    assert p["self_consistent_Hubble_from_dark_sector"] is False
    assert p["fixed_SM_entropy_degrees"]
