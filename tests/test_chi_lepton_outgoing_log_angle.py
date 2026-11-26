"""Finite-mediator CM angular mapping is exact for I0/I1, not full 5D."""
import pytest
from kinetics.chi_lepton_outgoing_log_angle import (
    exact_propagator_moments,cm_angle_moments,
)
from kinetics.independent_chi_lepton_5d_quadrature import evaluate


@pytest.mark.parametrize("A",(.01,.1,1,10,100,1000,10000))
def test_log_transfer_resolves_analytic_transport_kernel(A):
    exact=exact_propagator_moments(A)
    got=cm_angle_moments(A,8,"log_transfer")
    for x,y in zip(got,exact):
        assert abs(x/y-1)<2e-8


def test_uniform_cm_gauss_misses_peaked_vector_while_log_map_recovers():
    A=1000.
    exact=exact_propagator_moments(A)
    bad=cm_angle_moments(A,8,"uniform_cos")
    good=cm_angle_moments(A,8,"log_transfer")
    assert bad[0]/exact[0]<.1
    assert all(abs(x/y-1)<1e-8 for x,y in zip(good,exact))
    literal=cm_angle_moments(A,8,"propagator_compensated")
    assert abs(literal[0]/exact[0]-1)<1e-13
    assert abs(literal[1]/exact[1]-1)>.3


def test_full_pauli_narrow_vector_log_angle_controls_outgoing_bias():
    kw=dict(M=.1,m=0.,mv=.01,nr=4,nmu=6,nphi=12,
            radial_scheme="segmented_p",segments=8)
    lo=evaluate("smooth",nc=8,angular_scheme="log_transfer",**kw)
    hi=evaluate("smooth",nc=32,angular_scheme="log_transfer",**kw)
    bad=evaluate("smooth",nc=32,angular_scheme="uniform_cos",**kw)
    assert lo["Qchi"]>0 and lo["entropy"]>0
    assert abs(lo["Qchi"]/hi["Qchi"]-1)<3e-6
    assert abs(lo["entropy"]/hi["entropy"]-1)<1e-5
    assert abs(bad["Qchi"]/hi["Qchi"]-1)>.15
    assert hi["negative_entropy_events"]==0


def test_full_nonthermal_rate_needs_other_independent_grids():
    kw=dict(M=.1,m=0.,mv=.1,nr=4,nmu=6,nphi=12,
            radial_scheme="segmented_p",segments=8,
            angular_scheme="log_transfer")
    a=evaluate("oscillatory",nc=8,**kw)
    b=evaluate("oscillatory",nc=16,**kw)
    assert a["entropy"]>0 and b["entropy"]>0
    assert abs(a["Qchi"]-b["Qchi"])>0
