"""Finite-capacity SM photon-like bath and elastic fermionic heat transfer."""
from __future__ import annotations
import math
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.special import expit

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from sm_radiation_heat_feedback import (
    comoving_radiation_energy,comoving_radiation_entropy,
    radiation_elastic_exchange_step,demonstration)


def initial():
    q=np.geomspace(1e-4,30,96)
    f=np.clip(expit(-np.hypot(q,.1))*
              (1+.05*np.tanh(np.log1p(q))),0,1)
    return q,f


def test_radiation_baseline_and_exact_ideal_scaling():
    E=comoving_radiation_energy(1.,1.,55.)
    assert math.isclose(E,math.pi**2/30*55,rel_tol=1e-14)
    assert math.isclose(comoving_radiation_energy(2.,.5,55.),E,rel_tol=1e-14)
    assert math.isclose(
        comoving_radiation_entropy(2.,.5,55.),
        comoving_radiation_entropy(1.,1.,55.),rel_tol=1e-14)


def test_elastic_SM_heat_exchange_preserves_actual_total_energy_and_number():
    q,f=initial()
    r=radiation_elastic_exchange_step(q,f)
    assert abs(r["mother_number_change"])<1e-11
    assert abs(r["total_comoving_collision_energy_residual"])<1e-11
    assert math.isclose(r["total_energy_before"],r["total_energy_after"],rel_tol=1e-12)
    assert r["SM_radiation_comoving_energy_after"]>r["SM_radiation_comoving_energy_before"]
    assert r["mother_energy_after"]<r["mother_energy_before"]
    assert r["T_SM_after"]>r["T_SM_before"]
    assert r["total_entropy_change"]>0
    assert r["fermionic_KL_after"]<r["fermionic_KL_before"]


def test_pauli_bounds_and_large_relaxation_rate_exact_step():
    q,f=initial()
    r=radiation_elastic_exchange_step(q,f,gamma_el=1e5,dt=1e-2)
    assert 0<=r["Pauli_occupation_min_after"]<=1
    assert 0<=r["Pauli_occupation_max_after"]<=1
    assert abs(r["mother_number_change"])<1e-11
    assert abs(r["total_comoving_collision_energy_residual"])<1e-11
    assert r["total_entropy_change"]>0


def test_zero_elapsed_time_zero_energy_transfer_and_no_temperature_change():
    q,f=initial()
    for kwargs in ({"dt":0.},{"gamma_el":0.}):
        r=radiation_elastic_exchange_step(q,f,**kwargs)
        assert np.array_equal(r["mother_before"],r["mother_after"])
        assert r["T_SM_after"]==r["T_SM_before"]
        assert r["SM_radiation_comoving_energy_after"]==r["SM_radiation_comoving_energy_before"]
        assert r["total_entropy_change"]==0


def test_bath_temperature_feedback_sign_inverts_for_cold_parent():
    q=np.geomspace(1e-4,30,96)
    f=expit(-10*(q+.1))
    r=radiation_elastic_exchange_step(q,f,T_bath=1.,gamma_el=3.,dt=.3)
    assert r["mother_energy_after"]>r["mother_energy_before"]
    assert r["T_SM_after"]<r["T_SM_before"]
    assert r["total_entropy_change"]>0


@pytest.mark.parametrize("bad",[
    {"a":0},{"T_bath":0},{"gstar":0},{"gstar":-1},
    {"gamma_el":-1},{"dt":-1},{"mass":0},{"g1":0},
])
def test_unphysical_energy_reservoir_or_scattering_inputs_fail_closed(bad):
    q,f=initial()
    with pytest.raises(ValueError):
        radiation_elastic_exchange_step(q,f,**bad)


def test_physical_BEST_integration_and_CLASS_are_not_falsely_claimed():
    r=demonstration()
    assert r["finite_temperature_feedback_implemented"]
    assert r["time_integrator_first_order_bath_feedback"]
    assert r["is_one_fixed_a_collision_substep"]
    assert not r["microphysical_SM_scattering_matrix_element_derived"]
    assert not r["full_FRW_pressure_and_1to2_scatterings_integrated"]
    assert not r["physical_86eV_relic_or_CLASS_ready"]
