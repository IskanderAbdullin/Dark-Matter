"""Finite spacelike mediator in the independent five-dimensional Pauli heat.

Run: python -m pytest -q tests/test_finite_mediator_5d.py
These are conditional fixed-temperature checks, NOT a physical gamma_el.
"""
import math
import pytest
from kinetics.independent_massive_chie_5d import cm_heat


def q(mv, n=12, tc=.9):
    M, m, T = .1, 0., 1.
    pmax=math.sqrt((M+18*T)**2-M*M)
    return cm_heat(M,m,.001,T,tc,pmax=pmax,kmax=18*T,
                   nr=n,nang=n,nphi=n,mediator_mass=mv)


def test_contact_limit_and_finite_mediator_effect():
    reference=q(None)
    assert reference>0
    assert q(1e6)==pytest.approx(reference,rel=1e-9)
    assert abs(q(1.)-reference)>1e-3*abs(reference)


def test_common_temperature_detailed_balance_and_invalid_mass():
    assert abs(q(1.,tc=1.))<1e-17
    with pytest.raises(ValueError):
        q(0.)


def test_finite_mediator_refinement_and_heat_sign():
    # No assertion of convergence to a known continuum value at this grid.
    for n in (12,16):
        assert math.isfinite(q(1.,n=n))
        assert q(1.,n=n)>0
