"""Physical moving-particle energy-moment limits and signed recoil."""
import math
import sys
from pathlib import Path
import pytest
from scipy.integrate import quad
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from microphysics_portal import vector_contact_dsdt,sigma_transport_massless_e
from microphysics_recoil_heat import recoil_heating_rate_chi_at_rest
from microphysics_moving_heat import (
  moving_chi_unblocked_energy_moment,energy_kick_scattering_azimuth_average,
  incoming_s_and_moller_flux)


@pytest.mark.parametrize("M,p,k,mu",[
    (1.,0.,.1,.7),(1.,.1,.2,-.6),(.3,.8,.4,.1),(.3,1.5,.8,.9)])
def test_exact_cm_azimuth_averaged_energy_kick_by_independent_t_integration(M,p,k,mu):
    s,_=incoming_s_and_moller_flux(M,p,k,mu)
    Estart=math.hypot(M,p)
    gamma_Estar=(Estart+k)*(s+M*M)/(2*s)
    coeff_direct=gamma_Estar-Estart
    coeff=energy_kick_scattering_azimuth_average(M,p,k,mu)
    assert math.isclose(coeff,coeff_direct,rel_tol=1e-10,abs_tol=1e-13)
    pstar2=(s-M*M)**2/(4*s)
    tmin=-4*pstar2
    independent=quad(lambda t: vector_contact_dsdt(s,t,M,.17)*
                      coeff_direct*(-t)/(2*pstar2),
                      tmin,0,epsabs=0,epsrel=1e-11)[0]
    formula=coeff*sigma_transport_massless_e(s,M,.17)
    assert math.isclose(independent,formula,rel_tol=1e-10,abs_tol=1e-15)


def test_moving_unblocked_energy_reduces_to_real_rest_two_body_recoil():
    M=1.;T=.1;G=.2
    rest=recoil_heating_rate_chi_at_rest(M,T,G,blocking=False)
    ours=moving_chi_unblocked_energy_moment(M,0.,T,G)
    slow=moving_chi_unblocked_energy_moment(M,.0001,T,G)
    assert rest>0 and math.isclose(ours,rest,rel_tol=1e-13)
    assert abs(slow/rest-1)<1e-4


def test_hot_bath_heats_slow_chi_and_fast_chi_loses_energy():
    M=1.;T=.1;G=.2
    a=moving_chi_unblocked_energy_moment(M,.1,T,G)
    b=moving_chi_unblocked_energy_moment(M,1.,T,G)
    assert a>0 and b<0
    assert math.isclose(moving_chi_unblocked_energy_moment(M,1.,T,2*G)/b,
                        4,rel_tol=1e-10)


def test_zero_coupling_and_off_shell_domain():
    assert moving_chi_unblocked_energy_moment(1.,.3,.1,0.)==0
    with pytest.raises(ValueError):
        moving_chi_unblocked_energy_moment(1.,-.1,.1,1.)
    with pytest.raises(ValueError):
        energy_kick_scattering_azimuth_average(1.,.1,.2,1.1)
