"""Independent recoil-energy numerical controls and heavy-chi asymptote."""
import math
import sys
from pathlib import Path
from scipy.integrate import quad
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from microphysics_portal import sigma_transport_massless_e, vector_contact_dsdt
from microphysics_recoil_heat import (
    energy_gain_per_incident_energy,recoil_heating_rate_chi_at_rest,
    heavy_chi_recoil_asymptote)


def test_unblocked_recoil_identity_vs_transport_cross_section():
    M=1.;G=.15;T=.1
    for k in (.01,.1,1.):
        s=M*M+2*M*k
        pstar2=(s-M*M)**2/(4*s)
        expect=pstar2/M*sigma_transport_massless_e(s,M,G)
        got=energy_gain_per_incident_energy(k,M,G,T,blocking=False)
        assert math.isclose(got,expect,rel_tol=3e-13)


def test_independently_integrated_recoil_angular_moment():
    M=.5;k=.2;G=.17;T=.12
    s=M*M+2*M*k;tmin=-(s-M*M)**2/s
    def e_final(t):
        ef=k+t/(2*M)
        return 1-1/(math.exp(ef/T)+1)
    ref=quad(lambda t:vector_contact_dsdt(s,t,M,G)*
              (-t)/(2*M)*e_final(t),tmin,0,epsabs=0,epsrel=1e-12)[0]
    got=energy_gain_per_incident_energy(k,M,G,T)
    assert math.isclose(ref,got,rel_tol=2e-12)


def test_blocking_suppresses_positive_recoil_heating():
    M=1.;T=.3;G=.25
    full=recoil_heating_rate_chi_at_rest(M,T,G,blocking=True)
    unblocked=recoil_heating_rate_chi_at_rest(M,T,G,blocking=False)
    assert 0<full<unblocked
    assert math.isclose(recoil_heating_rate_chi_at_rest(M,T,2*G),4*full,
                        rel_tol=2e-12)


def test_heavy_rest_limit_with_pauli_blocking():
    M=100.;T=.001;G=.2
    real=recoil_heating_rate_chi_at_rest(M,T,G)
    asym=heavy_chi_recoil_asymptote(M,T,G)
    assert real>0 and abs(real/asym-1)<.001


def test_zero_coupling_and_unphysical_inputs():
    assert recoil_heating_rate_chi_at_rest(1.,.2,0.)==0
    assert energy_gain_per_incident_energy(.2,1.,0.,.2)==0
    with pytest.raises(ValueError):
        recoil_heating_rate_chi_at_rest(0.,.2,1.)
    with pytest.raises(ValueError):
        energy_gain_per_incident_energy(-.1,1.,1.,.2)
