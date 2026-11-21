"""Finite-q conservative elastic chi1 bath projection, quantum FD and energy Q."""
from __future__ import annotations
import math
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.special import expit
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from number_conserving_sm_bath import (
    number_conserving_bath_rates, exact_frozen_bath_step,demonstration)


def fixture():
    q=np.geomspace(1e-4,30.,96)
    initial=np.clip(expit(-np.hypot(q,.1))*
                    (1+.05*np.tanh(np.log1p(q))),0,1)
    return q,initial


def test_fd_chemical_mu_projects_elastic_number_not_externally_forces_mu_zero():
    q,_=fixture()
    a=1.18; mass=.1;T=.82; mu=.12
    f=expit((mu-np.hypot(q/a,mass))/T)
    result=exact_frozen_bath_step(q,f,a=a,mass=mass,Tbath=T)
    assert math.isclose(result["chemical_mu_number_projected"],mu,
                        rel_tol=1e-10,abs_tol=1e-10)
    assert np.max(abs(result["after"]-f))<1e-10
    assert abs(result["dark_comoving_energy_change"])<1e-10
    assert result["bath_T_not_self_consistently_evolved"]


def test_exact_bgk_step_preserves_mother_particle_number_and_pauli_bound():
    q,f=fixture()
    r=exact_frozen_bath_step(q,f,gamma=1.2,dt=.4)
    assert abs(r["comoving_number_residual"])<1e-11
    assert r["minimum_occupation"]>=0
    assert r["maximum_occupation"]<=1
    assert r["fermionic_KL_before"]>0
    assert r["fermionic_KL_after"]<r["fermionic_KL_before"]
    assert r["elastic_chi1_number_conserved"]


def test_sm_energy_sink_balances_actual_dark_transfer_with_frozen_a():
    q,f=fixture()
    r=exact_frozen_bath_step(q,f,gamma=2.,dt=.5,a=1.18)
    assert abs(r["dark_comoving_energy_change"])>1e-6
    assert math.isclose(
        r["SM_bath_comoving_energy_change"],
        -r["dark_comoving_energy_change"],rel_tol=1e-12)
    assert r["dark_plus_bath_comoving_energy_exchange_residual"]==0
    assert r["not_a_derived_SM_2to2_matrix_element"]


def test_bgk_instantaneous_gain_minus_loss_has_zero_number_source():
    q,f=fixture()
    rate=number_conserving_bath_rates(q,f,gamma=1.2,a=1.18)
    assert np.all(rate["gain"]>=0)
    assert np.all(rate["loss"]>=0)
    assert np.allclose(rate["rate"],rate["gain"]-rate["loss"])
    assert abs(rate["comoving_number_derivative"])<1e-10
    assert math.isclose(
        rate["comoving_dark_energy_derivative"],
        -rate["comoving_SM_bath_energy_derivative"],rel_tol=1e-12)


def test_very_stiff_exact_bgk_relaxation_still_positive_and_number_conserving():
    q,f=fixture()
    r=exact_frozen_bath_step(q,f,gamma=1e8,dt=1.)
    assert abs(r["comoving_number_residual"])<1e-11
    assert np.min(r["after"])>=0
    assert np.max(r["after"])<=1
    assert np.max(abs(r["after"]-r["target_fermi_dirac"]))<1e-12
    assert r["fermionic_KL_after"]<r["fermionic_KL_before"]


def test_zero_duration_or_zero_rate_preserves_physical_distribution():
    q,f=fixture()
    for kwargs in ({"dt":0},{"gamma":0}):
        r=exact_frozen_bath_step(q,f,**kwargs)
        assert np.array_equal(r["before"],r["after"])
        assert r["dark_comoving_energy_change"]==0


@pytest.mark.parametrize("bad",[
    {"a":0},{"mass":0},{"Tbath":0},{"g":0},
    {"gamma":-1},{"dt":-1},{"gamma":float("nan")}
])
def test_unphysical_bath_parameters_fail_closed(bad):
    q,f=fixture()
    with pytest.raises(ValueError):
        exact_frozen_bath_step(q,f,**bad)


def test_bad_distribution_and_momentum_mesh_fail_closed():
    q,f=fixture()
    for malformed_q,malformed_f in (
        (q[::-1],f),(q,f*3),(q,np.full_like(q,float("nan"))),
        (q[:3],f[:3]),(q,np.ones(len(q)-1))):
        with pytest.raises(ValueError):
            exact_frozen_bath_step(malformed_q,malformed_f)


def test_demonstration_does_not_fake_plasma_temperature_feedback_or_dm_relic():
    r=demonstration()
    assert r["fermionic_KL_after"]<r["fermionic_KL_before"]
    assert r["bath_T_not_self_consistently_evolved"]
    assert r["not_a_converged_cosmological_DM_model"]
    assert abs(r["comoving_number_residual"])<1e-11
