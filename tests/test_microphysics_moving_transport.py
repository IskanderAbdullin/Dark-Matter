"""Independent limiting cases of moving-chi incoming collision frequency."""
import math
import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from microphysics_moving_transport import moving_chi_transport_frequency
from microphysics_portal import bath_transport_frequency_chi_at_rest


def test_exact_rest_reduction_and_small_nonzero_momentum():
    M=.2;T=.1;G=.6
    ref=bath_transport_frequency_chi_at_rest(M,T,G)
    assert math.isclose(moving_chi_transport_frequency(M,0.,T,G),ref,rel_tol=1e-13)
    almost=moving_chi_transport_frequency(M,M*1e-4,T,G)
    assert abs(almost/ref-1)<1e-6


def test_moving_frequency_positive_finite_and_coupling_quadratic():
    M=.2;T=.1;G=.3
    for p in (.01,.2,2.):
        val=moving_chi_transport_frequency(M,p,T,G)
        assert val>0 and math.isfinite(val)
        assert math.isclose(moving_chi_transport_frequency(M,p,T,2*G)/val,4,
                            rel_tol=1e-11)


def test_no_coupling_and_invalid_domain():
    assert moving_chi_transport_frequency(.2,.5,.1,0.)==0
    for args in [(.2,-.01,.1,.3),(.2,.1,0.,.3),
                 (.2,.1,.1,-.1),(.2,float("nan"),.1,.3)]:
        with pytest.raises(ValueError):
            moving_chi_transport_frequency(*args)
