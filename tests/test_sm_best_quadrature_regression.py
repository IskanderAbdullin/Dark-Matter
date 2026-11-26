"""Critical BEST Simpson-in-log-q vs elastic SM BGK trapezoid mismatch.

The first *real* coupled BEST+SM experiment conserved a BGK's own
trapezoidal number but drifted in BEST's published log-grid Simpson
moments. The BH source is imported only in GitHub Actions; reproduce
its documented quadrature algorithm here for independent CI tests.
"""
from __future__ import annotations
import math
import sys
from pathlib import Path
import numpy as np
import pytest
from scipy.special import expit

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from number_conserving_sm_bath import (
    exact_frozen_bath_step,number_conserving_bath_rates)
from sm_radiation_heat_feedback import radiation_elastic_exchange_step


def best_log_simpson_dq(q):
    """Independent log-Simpson/3-8 weights matching pinned BEST v1.3.5."""
    q=np.asarray(q,float)
    intervals=len(q)-1
    h=np.diff(np.log(q)).mean()
    m=intervals if intervals%2==0 else intervals-3
    w=np.zeros_like(q)
    w[:m+1]=h/3
    w[1:m:2]*=4
    w[2:m:2]*=2
    if m<intervals:
        w[m:]+=3*h/8*np.array([1.,3.,3.,1.])
    return w*q


def testcase():
    q=np.geomspace(1e-4,30.,56)
    f=np.clip(expit(-np.hypot(q,.1))*
              (1+.05*np.tanh(np.log1p(q))),0,1)
    return q,f,best_log_simpson_dq(q)


def best_n_and_E(q,f,dq,a=1.18,mass=.1,g=4):
    W=g/(2*np.pi**2)*dq*q**2
    return (float(np.dot(W,f)),float(np.dot(W*np.hypot(q,a*mass),f)))


def test_old_trapezoid_sm_contact_did_not_preserve_actual_best_number_moment():
    q,f,dq=testcase()
    old=exact_frozen_bath_step(q,f,a=1.18,mass=.1,
                               Tbath=.82,gamma=100000,dt=2.5e-6)
    n_old=best_n_and_E(q,f,dq)[0]
    n_after=best_n_and_E(q,old["after"],dq)[0]
    assert abs(old["comoving_number_residual"])<1e-11
    # Regression witness of the REAL historical simulation defect.
    assert abs(n_after-n_old)>1e-8


def test_injected_best_weights_preserve_number_and_full_chi1_energy_plus_SM():
    q,f,dq=testcase()
    r=radiation_elastic_exchange_step(
        q,f,a=1.18,mass=.1,T_bath=.82,g1=4,
        gamma_el=100000,dt=2.5e-6,gstar=55,
        quadrature_weights=dq)
    assert r["quadrature_matched_to_supplied_external_solver"]
    n0,e0=best_n_and_E(q,r["mother_before"],dq)
    n1,e1=best_n_and_E(q,r["mother_after"],dq)
    assert abs(n1-n0)<1e-11
    assert abs(r["mother_number_change"])<1e-11
    assert math.isclose(r["mother_energy_before"],e0,rel_tol=1e-14)
    assert math.isclose(r["mother_energy_after"],e1,rel_tol=1e-14)
    assert math.isclose(r["SM_radiation_comoving_energy_after"]-
                        r["SM_radiation_comoving_energy_before"],
                        e0-e1,rel_tol=1e-8,abs_tol=1e-12)
    assert abs(r["total_comoving_collision_energy_residual"])<1e-11
    assert r["total_entropy_change"]>=-1e-12
    assert 0<=r["Pauli_occupation_min_after"]<=r["Pauli_occupation_max_after"]<=1


def test_injected_best_weighted_rate_has_zero_exact_discrete_number_source():
    q,f,dq=testcase()
    x=number_conserving_bath_rates(
        q,f,a=1.18,mass=.1,Tbath=.82,gamma=100000,
        quadrature_weights=dq)
    assert abs(x["comoving_number_derivative"])<1e-6
    assert math.isclose(x["comoving_dark_energy_derivative"],
                        -x["comoving_SM_bath_energy_derivative"],
                        rel_tol=1e-13)


def test_invalid_weight_shapes_negatives_and_nan_rejected_without_fallback():
    q,f,dq=testcase()
    for bad in (dq[:-1],np.full_like(dq,np.nan),
                np.zeros_like(dq),-dq):
        with pytest.raises(ValueError,match="quadrature_weights"):
            exact_frozen_bath_step(q,f,quadrature_weights=bad)
        with pytest.raises(ValueError,match="quadrature_weights"):
            radiation_elastic_exchange_step(q,f,quadrature_weights=bad)


def test_new_method_uses_explicit_external_weights_not_trapezoid_silently():
    q,f,dq=testcase()
    x=exact_frozen_bath_step(q,f,quadrature_weights=dq)
    y=exact_frozen_bath_step(q,f)
    assert abs(x["comoving_number_residual"])<1e-11
    assert abs(y["comoving_number_residual"])<1e-11
    assert abs(x["chemical_mu_number_projected"]-
               y["chemical_mu_number_projected"])>1e-7
