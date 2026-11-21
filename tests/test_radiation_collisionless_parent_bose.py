"""Independent radiation-Bose parent-redshift analytic benchmark and prior BEST."""
from __future__ import annotations
import json
import math
import sys
from pathlib import Path

import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from radiation_collisionless_parent_bose import (
    scale_factor,collisionless_parent_fd_at_physical_energy,
    gain_coefficient,frozen_parent_rd_occupation,
    frozen_parent_rd_ir_number,evaluate)
from massive_daughter_decay_source import zero_momentum_bose_occupation_source
from best_early_bose_exponential import initial_gain_coefficient


def test_rd_background_and_comoving_fd_mother_not_thermally_reset():
    assert math.isclose(scale_factor(.00014,.0001),math.sqrt(1.4),rel_tol=1e-14)
    # At a=1 the exact occupation equals FD of the physical mother energy.
    E=5.
    assert math.isclose(
        collisionless_parent_fd_at_physical_energy(E,1.,.1),
        1/(1+math.exp(E)),rel_tol=1e-13)
    # At a>1 the same physical energy is rarer: its COMOVING momentum
    # is larger, mother is not reset to thermal FD at its current mass.
    assert collisionless_parent_fd_at_physical_energy(E,1.2,.1)<\
           collisionless_parent_fd_at_physical_energy(E,1.,.1)


@pytest.mark.parametrize("q",[0.,1e-6,2e-5,1e-4])
def test_initial_mother_gain_agrees_with_independent_fixed_a_source(q):
    a=gain_coefficient(q,t=1e-4)
    b=initial_gain_coefficient(q)
    assert math.isclose(a,b,rel_tol=1e-11,abs_tol=1e-9)


def test_exact_massive_zero_momentum_initial_gain_limit():
    q0=gain_coefficient(0.,1e-4)
    assert math.isclose(q0,zero_momentum_bose_occupation_source(),
                        rel_tol=1e-12)
    assert 26000<q0<28000


def test_actual_rd_mother_redshift_suppresses_soft_gain():
    q=1e-4
    begin=gain_coefficient(q,1e-4)
    end=gain_coefficient(q,1.4e-4)
    assert math.isclose(begin,27046.973477458418,rel_tol=1e-10)
    assert math.isclose(end,10928.943290653331,rel_tol=1e-10)
    assert end<begin/2


def test_short_time_time_integrated_bose_exponential_matches_reference():
    result=evaluate()
    assert math.isclose(result["fphi_at_qmin"],1.0196839876614525,
                        rel_tol=2e-9)
    assert math.isclose(result["fphi_at_qzero"],1.000859744389198,
                        rel_tol=2e-9)
    assert math.isclose(result["below_qmin_comoving_number"],
                        1.7092584190101826e-14,rel_tol=2e-7)
    assert not result["macroscopic_condensate_or_real_86eV_FRW_relic_or_CLASS_validated"]


def test_analytic_rd_explanation_agrees_with_preliminary_best_time_extrapolation():
    path=(Path(__file__).resolve().parents[1]/
          "kinetics/results/best-rd-time-refinement-2026-11-21.json")
    data=json.loads(path.read_text(encoding="utf-8"))
    best_first_order_extrapolated=data["interpreted_convergence_order"][
        "preliminary_Richardson_f_qmin_16_32"]
    analytic=frozen_parent_rd_occupation(1e-4)
    assert abs(best_first_order_extrapolated/analytic-1)<.002
    # This is NOT proof of simultaneous time/q/MC convergence.
    assert data["claimed_time_convergence"] is False


@pytest.mark.parametrize("kwargs",[
    {"t":.00009},{"q":-.01},{"mphi":0.},{"M":0.},
    {"gphi":0.},{"rest_width":0.},{"q":float("nan")}
])
def test_invalid_cosmic_time_and_kinematics_fail_closed(kwargs):
    with pytest.raises(ValueError):
        gain_coefficient(**{"q":.0001,"t":.00011,**kwargs})


def test_invalid_cutoff_and_negative_seed_rejected():
    with pytest.raises(ValueError):
        frozen_parent_rd_ir_number(q_cut=0)
    with pytest.raises(ValueError):
        frozen_parent_rd_occupation(.001,seed=-1)
