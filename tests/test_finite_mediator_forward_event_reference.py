"""Independent finite-vector mediator polar angular reference and negative grid control."""
import math
from kinetics.finite_mediator_forward_event_reference import (
    angular_integrals,propagator_angular_identity,
)


def test_closed_propagator_jacobian_at_many_forward_scales():
    for ps2 in (.001,.1,2.):
        for mV in (.003,.03,.3,3.,100.):
            a,b=propagator_angular_identity(ps2,mV)
            assert math.isclose(a,b,rel_tol=5e-15)


def test_fixed_event_light_mediator_not_converged_at_48_nodes():
    d=angular_integrals(.03,orders=(24,48,96,160))
    vals={r["order"]:r for r in d["scan"]}
    assert d["quad_reference"]>0
    assert abs(vals[48]["cosine_relative_error"])>.01
    assert abs(vals[48]["transfer_relative_error"])>.005
    assert abs(vals[160]["transfer_relative_error"])<4e-5
    assert abs(vals[160]["cosine_relative_error"])<7e-5


def test_fixed_event_moderate_mediator_agrees_with_adaptive_quad():
    d=angular_integrals(.3,orders=(16,24,48))
    vals={r["order"]:r for r in d["scan"]}
    assert abs(vals[16]["transfer_relative_error"])<1e-4
    assert abs(vals[48]["cosine_relative_error"])<1e-10
    assert abs(vals[48]["transfer_relative_error"])<1e-10
