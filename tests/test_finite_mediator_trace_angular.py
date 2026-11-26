"""Independent numerical quadrature of nonconstant finite-vector spin trace."""
import math
import pytest
from scipy.integrate import quad
from kinetics.finite_mediator_trace_angular import moments, trace_angular


@pytest.mark.parametrize("b",[0.,1e-12,1e-5,.01,.1,1.,10.,100.])
def test_moments_against_adaptive_quadrature(b):
    for n, val in enumerate(moments(b)):
        ref=quad(lambda y:y**n/(1+b*y)**2,0.,2.,epsabs=1e-13)[0]
        assert val == pytest.approx(ref,rel=3e-11,abs=1e-13)


@pytest.mark.parametrize("M,m,ps,mV",[
    (.1,0.,.01,.02),(.1,0.,1.,.2),(.1,0.,1.,10.),
    (.1,.00051099895,.3,.01),(.1,.00051099895,.3,100.),
    (1.,.5,0.,2.),(1.,.5,1.,1.)
])
def test_full_spin_trace_against_direct_t_integral(M,m,ps,mV):
    e1=math.hypot(M,ps)
    e2=math.hypot(m,ps)
    s=(e1+e2)**2
    D=s-M*M-m*m
    B=M*M+m*m
    a=2*ps*ps
    ref=quad(lambda y:(2*D*D-2*a*(D+B)*y+a*a*y*y)/
             (1+a*y/mV**2)**2,0.,2.,epsabs=1e-12)[0]
    assert trace_angular(s,M,m,mV)==pytest.approx(ref,rel=3e-11,abs=1e-12)
