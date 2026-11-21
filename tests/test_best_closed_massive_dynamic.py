"""A future actual BEST 3-species kinetic evolution cannot be faked by protocol."""
from __future__ import annotations
import math
import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from best_closed_massive_dynamic import dynamic_protocol,execute_best
from best_external_vacuum_reference import BEST_PIN


def test_default_dynamic_bose_bath_protocol_has_no_external_run():
    p=dynamic_protocol()
    assert p["massive_initial_vacuum_protocol"]["verified_best_git_commit_expected"]==BEST_PIN
    assert p["steps"]==5
    assert p["dtGamma_per_step"]==1e-5
    assert math.isclose(p["total_tGamma"],5e-5,rel_tol=1e-15)
    assert p["integrator"]=="exprb"
    assert "NOT EXECUTED" in p["status"]
    assert p["not_a_physical_86eV_FRW_or_CLASS"]
    assert p["daughter_initial_occupation"]==1e-12
    assert p["daughter_initial_seed_shape"].startswith("1e-12 * exp(-q)")
    assert p["static_scale_factor"]==1.


@pytest.mark.parametrize("opts",[
    {"steps":0},{"dtGamma":0},{"method":"unsupported"},
    {"mphi":.2},{"nq":1},{"vegas_neval":200}
])
def test_unphysical_or_unsupported_dynamic_protocol_rejected(opts):
    with pytest.raises(ValueError):
        dynamic_protocol(**opts)


def test_missing_best_source_cannot_produce_fake_time_series(tmp_path):
    p=dynamic_protocol(steps=2)
    with pytest.raises(ValueError,match="missing"):
        execute_best(tmp_path,tmp_path/"fake.json",p)
    assert not (tmp_path/"fake.json").exists()


def test_dynamic_protocol_has_physical_massive_zero_mode_context_but_no_BEC_claim():
    p=dynamic_protocol()
    assert p["massive_initial_vacuum_protocol"]["mass_phi"]>0
    assert "distinguish" in p["zero_mode"]
    assert p["massive_initial_vacuum_protocol"]["qmin"]==1e-4
    assert dynamic_protocol(qmin=5e-5)["massive_initial_vacuum_protocol"]["qmin"]==5e-5
    with pytest.raises(ValueError):
        dynamic_protocol(qmin=0)
    assert p["real_BEST_physics"].startswith("time evolving")
