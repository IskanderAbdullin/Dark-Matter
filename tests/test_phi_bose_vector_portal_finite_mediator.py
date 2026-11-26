"""Gauge-consistent χ1 AND χ2 finite-mass vector portal candidate to SM e±.

Differentiate elastic event transport frequency proxy vs actual thermal drag.
Tests include independent 4x4 Dirac gamma trace, known contact analytic
reference, massless MB thermal rate integral, finite mediator suppression,
physical q1=q2 Yukawa real-scalar gauge-charge requirement.
"""
import math
import numpy as np
from scipy.integrate import quad
from scipy.special import kve

from kinetics.microphysics_portal import (
    sigma_transport_massless_e,vector_contact_m2,
)
from kinetics.gauge_charge_yukawa_audit import analyze as gauge
from kinetics.phi_bose_vector_portal_finite_mediator import (
    ALPHA_EM, vector_portal_setup, amplitude_squared,
    finite_vector_cross_section, MB_transport_collision_proxy,
    ultrarelativistic_MB_contact_limit, t_bounds_massless_e,
)


def _dirac_matrices():
    I=np.eye(2,dtype=complex)
    Z=np.zeros((2,2),complex)
    sig=[np.array([[0,1],[1,0]],complex),
         np.array([[0,-1j],[1j,0]],complex),
         np.diag([1,-1]).astype(complex)]
    return [np.block([[I,Z],[Z,-I]])]+[
        np.block([[Z,a],[-a,Z]]) for a in sig]


def _gamma_independent_amp2(s,c,mchi=.02,gD=.1,epsilon=1e-5,mV=10.):
    """4x4 COMPLEX gamma trace with correctly lowered electron current."""
    rt=math.sqrt(s)
    Ep=(s+mchi*mchi)/(2*rt)
    pe=(s-mchi*mchi)/(2*rt)
    px=pe*math.sqrt(1-c*c)
    P=np.array([Ep,0.,0.,pe])
    K=np.array([pe,0.,0.,-pe])
    Pp=np.array([Ep,px,0.,pe*c])
    Kp=np.array([pe,-px,0.,-pe*c])
    t=-2*pe*pe*(1-c)
    gam=_dirac_matrices()
    lower=[(1 if mu==0 else -1)*g for mu,g in enumerate(gam)]
    def slash(v):
        return sum((sgn*gam[i]*v[i] for i,sgn in enumerate((1,-1,-1,-1))),
                   start=np.zeros((4,4),complex))
    I=np.eye(4,dtype=complex)
    T=0j
    for mu in range(4):
        for nu in range(4):
            chi=np.trace((slash(Pp)+mchi*I)@gam[mu]@
                          (slash(P)+mchi*I)@gam[nu])
            ele=np.trace(slash(Kp)@lower[mu]@slash(K)@lower[nu])
            T+=chi*ele
    product=gD*epsilon*math.sqrt(4*math.pi*ALPHA_EM)
    return float((T/4*product*product/(mV*mV-t)**2).real),t


def test_real_phi_Yukawa_requires_equal_chi1_chi2_vector_portal_charges():
    yes=gauge(1,1,0,real_phi=True)
    no=gauge(1,0,0,real_phi=True)
    assert yes["renormalizable_yukawa_gauge_invariant"]
    assert yes["chi1_diagonal_vector_current_nonzero"]
    assert yes["chi2_diagonal_vector_current_nonzero"]
    assert not no["renormalizable_yukawa_gauge_invariant"]
    assert vector_portal_setup()["q_chi1"]==vector_portal_setup()["q_chi2"]


def test_full_finite_mass_vector_photon_spin_average_against_independent_dirac_trace():
    for mchi in (.02,.1):
        for E in (.2,.5,2.):
            for c in (-.7,0.,.6):
                s=E*E
                if s<=mchi*mchi:continue
                exact,t=_gamma_independent_amp2(s,c,mchi)
                amp=amplitude_squared(s,t,mchi)
                assert exact>0 and amp>0
                assert math.isclose(amp,exact,rel_tol=4e-13),(mchi,E,c,exact,amp)


def test_full_t_channel_transport_contact_limit_and_finite_mediator_effect():
    s=25.
    m=.02
    setup=vector_portal_setup()
    c=sigma_transport_massless_e(s,m,setup["G_contact_GeV_minus2"])
    finite=finite_vector_cross_section(s,m)
    assert math.isclose(finite/c,.7535901224957986,rel_tol=3e-8)
    assert 0<finite<c
    for m in (.02,.1):
        for s in (.25,1.,25.):
            g=vector_portal_setup(mV=1000.)["G_contact_GeV_minus2"]
            orig=sigma_transport_massless_e(s,m,g)
            new=finite_vector_cross_section(s,m,mV=1000.)
            assert 0<new<=orig
            assert math.isclose(new/orig,1.,rel_tol=1e-4)


def test_MB_transport_collision_proxy_and_absolute_Hubble_1GeV_benchmark():
    # This is NOT actual kinetic drag or χ+SM plasma equilibration.
    r2=MB_transport_collision_proxy(1.,mchi=.02)
    r1=MB_transport_collision_proxy(1.,mchi=.1)
    assert math.isclose(r2["MB_transport_proxy_over_H_NOT_drag"],
                        2.5304567019732462,rel_tol=3e-5)
    assert math.isclose(r1["MB_transport_proxy_over_H_NOT_drag"],
                        2.5310687667878193,rel_tol=3e-5)
    for T,ref in ((.1,.0038588321694743646),
                  (.5,.4237213977788347),
                  (2.,10.43040273922921)):
        a=MB_transport_collision_proxy(T,mchi=.02)
        assert math.isclose(a["MB_transport_proxy_over_H_NOT_drag"],ref,
                            rel_tol=4e-5)


def test_MB_reaction_density_independent_sqrt_s_integration_normalization():
    T=1.
    m=.02
    # Different independent variable and Bessel K1 from u=sqrt(s)/T:
    # γ=gχ ge T/(32π4) int 2 dE σ_T(E²)(E²−m²)² K1(E/T).
    integrand=lambda E:(
        2*finite_vector_cross_section(E*E,m)*
        (E*E-m*m)**2*kve(1,E/T)*math.exp(-E/T))
    integral=quad(integrand,m,40*T,epsabs=0,epsrel=2e-5,limit=300)[0]
    gamma=(2*4)*T/(32*math.pi**4)*integral
    r=MB_transport_collision_proxy(T,m)
    assert math.isclose(gamma,
                        r["MB_transport_weighted_forward_scattering_density_GeV4"],
                        rel_tol=5e-5)
    # Repeat the independent sqrt(s) variable with massive χ1, not just χ2.
    m1=.1
    integral1=quad(lambda E:2*finite_vector_cross_section(E*E,m1)*
                   (E*E-m1*m1)**2*kve(1,E/T)*math.exp(-E/T),
                   m1,40*T,epsabs=0,epsrel=2e-5,limit=300)[0]
    gamma1=8*T/(32*math.pi**4)*integral1
    ref1=MB_transport_collision_proxy(T,mchi=m1)
    assert math.isclose(gamma1,
                        ref1["MB_transport_weighted_forward_scattering_density_GeV4"],
                        rel_tol=5e-5)



def test_analytic_ultrarelativistic_MB_contact_limit_not_confused_with_drag():
    m=.0001;T=1.;mv=1000.
    r=MB_transport_collision_proxy(T,mchi=m,mV=mv)
    ref=ultrarelativistic_MB_contact_limit(
        T,r["G_contact_GeV_minus2"])
    assert r["MB_transport_weighted_per_chi_collision_frequency_GeV"]>0
    assert math.isclose(
        r["MB_transport_weighted_per_chi_collision_frequency_GeV"]/ref,
        1.,rel_tol=2e-3)
