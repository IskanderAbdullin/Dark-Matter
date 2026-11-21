"""Early frozen-parent exponential is approximate but independently predictive."""
from __future__ import annotations
import math
import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from best_early_bose_exponential import (
    initial_gain_coefficient,frozen_bath_phi_occupation,
    frozen_bath_ir_number,evaluate
)
from massive_daughter_decay_source import zero_momentum_bose_occupation_source


def test_qzero_exact_gain_matches_separately_derived_massive_two_body_limit():
    zero=initial_gain_coefficient(0.)
    assert math.isclose(zero,zero_momentum_bose_occupation_source(),rel_tol=1e-14)
    assert 26000<zero<28000


def test_frozen_bath_local_bose_gain_exponent_exact():
    q=.0001;t=4e-5
    A=initial_gain_coefficient(q)
    f=frozen_bath_phi_occupation(q,t,seed=0)
    assert math.isclose(f,math.expm1(A*t),rel_tol=1e-14)
    assert math.isclose(f,1.9502176258293882,rel_tol=1e-9)
    assert f>A*t


def test_frozen_bath_predicted_below_qmin_particle_number_includes_bose_stimulation():
    t=4e-5
    full=frozen_bath_ir_number(qcut=1e-4,t=t,seed=0)
    assert math.isclose(full,3.2705325738051216e-14,rel_tol=2e-8)
    assert full>1.8e-14
    assert full<4e-14


def test_low_q_early_growth_is_linear_and_qvolume_exact_at_small_times():
    for q in (0.,1e-4,.001):
        t=1e-8
        A=initial_gain_coefficient(q)
        f=frozen_bath_phi_occupation(q,t,seed=0)
        assert abs(f/(A*t)-1)<.0002


def test_bose_exponential_cold_f_peak_is_not_equal_to_total_relic_number():
    r=evaluate()
    assert r["nonlinear_frozen_bath_phi_f_qmin"]>1
    assert r["nonlinear_frozen_bath_phi_q_below_qmin_number"]<1e-12
    assert r["has_no_explicit_condensate_equation"]
    assert r["not_self_consistent_FRW_or_CLASS"]


@pytest.mark.parametrize("bad",[
    {"t":-1},{"seed":-1},{"t":float("nan")}
])
def test_invalid_short_time_inputs_rejected(bad):
    with pytest.raises(ValueError):
        frozen_bath_phi_occupation(.001,**{"t":1e-5,**bad})


def test_invalid_negative_q_cutoff_is_rejected():
    with pytest.raises(ValueError):
        frozen_bath_ir_number(qcut=0)
