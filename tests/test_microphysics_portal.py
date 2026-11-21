"""Independent numerical and analytic checks for the deliberately limited EFT pilot."""
import math
import sys
from pathlib import Path
import pytest
from scipy.integrate import quad
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from microphysics_portal import (
  kallen_sqrt_decay,scalar_yukawa_spin_averaged_m2,scalar_yukawa_width,
  best_spin_summed_m2,width_from_best_spin_sum,vector_contact_m2,
  vector_contact_dsdt,sigma_total_massless_e,sigma_transport_massless_e,
  bath_transport_frequency_chi_at_rest,bath_transport_frequency_heavy_limit,
  radiation_H,exploratory_snapshot)


def test_massless_yukawa_trace_and_width_independent_two_body_formula():
    M=.1;y=2e-5
    expected=y*y*M/(16*math.pi)
    assert math.isclose(scalar_yukawa_width(M,0,0,y),expected,rel_tol=1e-14)
    amp=scalar_yukawa_spin_averaged_m2(M,0,0,y)
    pstar=M/2
    assert math.isclose(scalar_yukawa_width(M,0,0,y),amp*pstar/(8*math.pi*M*M),
                        rel_tol=1e-14)


@pytest.mark.parametrize("m2,mphi",[(0,.000001),(.007,.014),(.02,.051)])
def test_best_normalization_equal_width_derived_from_tree_lagrangian(m2,mphi):
    M=.1;y=1.1e-4;g1=4
    width=scalar_yukawa_width(M,m2,mphi,y)
    amplitude=best_spin_summed_m2(M,m2,mphi,y,g1)
    pstar=kallen_sqrt_decay(M,m2,mphi)/(2*M)
    assert math.isclose(amplitude,g1*8*math.pi*M*M*width/pstar,rel_tol=2e-14)
    assert math.isclose(width_from_best_spin_sum(M,m2,mphi,amplitude,g1),
                        width,rel_tol=2e-14)


def test_threshold_not_fixed_width_and_zero_coupling():
    assert scalar_yukawa_width(.1,.06,.04,1e-2)==0
    assert scalar_yukawa_width(.1,.06,.05,1e-2)==0
    assert scalar_yukawa_width(.1,.01,.02,0)==0
    assert scalar_yukawa_width(.1,.01,.02,1e-2)>0
    with pytest.raises(ValueError):
        scalar_yukawa_width(.1,-.01,.02,1e-2)
    with pytest.raises(ValueError):
        scalar_yukawa_width(.1,0,0,float("nan"))


@pytest.mark.parametrize("s,m,G",[(1.21,1.,.01),(.025,.1,1.2),(5.,.01,.2)])
def test_transport_and_total_integrals_against_independent_numerical_t_quadrature(s,m,G):
    tmin=-(s-m*m)**2/s
    p2=(s-m*m)**2/(4*s)
    cross=quad(lambda t:vector_contact_dsdt(s,t,m,G),tmin,0,
               epsabs=1e-12,epsrel=1e-12)[0]
    tran=quad(lambda t:(-t)/(2*p2)*vector_contact_dsdt(s,t,m,G),tmin,0,
              epsabs=1e-12,epsrel=1e-12)[0]
    assert math.isclose(cross,sigma_total_massless_e(s,m,G),rel_tol=4e-13)
    assert math.isclose(tran,sigma_transport_massless_e(s,m,G),rel_tol=4e-13)
    assert vector_contact_m2(s,0,m,G)>=0
    assert vector_contact_m2(s,tmin,m,G)>=0


def test_transport_low_and_high_energy_angular_limits():
    M=1.
    for eps,ratio in [(1e-5,2/3),(1e5,7/8)]:
        s=M*M*(1+eps)
        got=sigma_transport_massless_e(s,M,.2)/sigma_total_massless_e(s,M,.2)
        assert abs(got-ratio)<1e-4


def test_rest_frame_bath_rate_exact_coupling_scaling_and_heavy_m_limit():
    T=.001;M=100.;G=0.2
    a=bath_transport_frequency_chi_at_rest(M,T,G)
    b=bath_transport_frequency_chi_at_rest(M,T,2*G)
    assert a>0 and math.isclose(b/a,4,rel_tol=1e-12)
    limit=bath_transport_frequency_heavy_limit(T,G)
    assert abs(a/limit-1)<2e-3
    assert bath_transport_frequency_chi_at_rest(M,T,0)==0


def test_finite_temperature_transport_frequency_is_not_bgk_gamma():
    snap=exploratory_snapshot()
    assert snap["no_thermal_Gamma_or_bgk_calibration"]
    assert snap["no_frozen_relic_or_BEC_or_CLASS"]
    assert snap["tree_vacuum_Gamma_rest"]>0
    assert snap["radiation_H_approx"]>0
    assert snap["rest_transport_frequency_over_H_NOT_gamma_el_over_H"]>0


@pytest.mark.parametrize("fn,args",[
   (sigma_total_massless_e, (1.,1.,1.)),
   (sigma_transport_massless_e,(0.,1.,1.)),
   (bath_transport_frequency_chi_at_rest,(1.,-1.,1.)),
   (radiation_H,(1.,0.)),
])
def test_fail_closed_bad_inputs(fn,args):
    with pytest.raises(ValueError):
        fn(*args)
