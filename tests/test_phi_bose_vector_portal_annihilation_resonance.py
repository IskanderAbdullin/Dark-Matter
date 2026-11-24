"""Conditional dark vector χχbar↔e+e− resonance: independent Dirac trace,
Breit-Wigner theta quadrature, MB detailed balance and heat.

Not a full FD/Pauli thermal plasma or observed scalar dark relic.
"""
import math
import numpy as np
from scipy.integrate import quad
from scipy.special import kve

from kinetics.phi_bose_vector_portal_finite_mediator import ALPHA_EM
from kinetics.phi_bose_vector_portal_annihilation_resonance import (
    MCHI1,MCHI2,vector_to_Dirac_width,minimal_vector_width,
    annihilation_sigma,inverse_ee_sigma,s_channel_thermodynamic_MB,
    net_MB_two_temperature_ledger,
)


def _gamma_matrices():
    i=np.eye(2,dtype=complex);z=np.zeros((2,2),complex)
    p=[np.array([[0,1],[1,0]],complex),
       np.array([[0,-1j],[1j,0]],complex),
       np.diag([1,-1]).astype(complex)]
    return [np.block([[i,z],[z,-i]])]+[
        np.block([[z,s],[-s,z]]) for s in p]


def _independent_gamma_trace_sigma(s,m,ctheta,mV=10.,gD=.1,epsilon=1e-5):
    E=math.sqrt(s)/2
    pin=math.sqrt(s/4-m*m)
    pout=E
    pc=np.array([E,0.,0.,pin])
    pa=np.array([E,0.,0.,-pin])
    sine=math.sqrt(1-ctheta*ctheta)
    ke=np.array([E,pout*sine,0.,pout*ctheta])
    ka=np.array([E,-pout*sine,0.,-pout*ctheta])
    G=_gamma_matrices()
    sig=[1,-1,-1,-1]
    def slash(p):
        return sum((sig[i]*g*p[i] for i,g in enumerate(G)),
                   start=np.zeros((4,4),complex))
    I=np.eye(4,dtype=complex)
    tr=0j
    for mu in range(4):
        for nu in range(4):
            ferm=np.trace((slash(pa)-m*I)@G[mu]@(slash(pc)+m*I)@G[nu])
            electron=np.trace(slash(ke)@(sig[mu]*G[mu])@slash(ka)@
                              (sig[nu]*G[nu]))
            tr+=ferm*electron
    gp2=(gD*epsilon*math.sqrt(4*math.pi*ALPHA_EM))**2
    width=minimal_vector_width(mV,gD,epsilon)["width_GeV"]
    D=(s-mV*mV)**2+(mV*width)**2
    amp2=gp2/D*float(tr.real/4)
    return amp2/(32*math.pi*s)*(pout/pin)


def test_exact_vacuum_s_channel_spin_average_vs_independent_4x4_Dirac_matrices():
    from numpy.polynomial.legendre import leggauss
    z,w=leggauss(16)
    for m in (MCHI1,MCHI2):
        for root_s in (.5,1.,3.,10.,15.):
            s=root_s**2
            reference=annihilation_sigma(s,m)
            independent=sum(ww*_independent_gamma_trace_sigma(
                s,m,float(zz)) for zz,ww in zip(z,w))
            assert reference>0
            assert math.isclose(reference,independent,rel_tol=2e-12),(
                m,root_s,reference,independent)


def test_finite_vector_width_is_fixed_by_same_both_Dirac_charges():
    w=minimal_vector_width(10.,.1,1e-5)
    assert math.isclose(w["width_GeV"],
                        .005305164612731479,rel_tol=3e-12)
    assert set(w["channels_GeV"])=={
        "chi1+antichi1","chi2+antichi2","e-+e+"}
    assert w["channels_GeV"]["e-+e+"]<1e-10
    assert w["channels_GeV"]["chi1+antichi1"]>0
    assert w["channels_GeV"]["chi2+antichi2"]>0
    assert vector_to_Dirac_width(.1,.1,.1)==0


def test_reversible_inverse_ee_cross_section_has_exact_Microreversible_flux():
    for m in (.02,.1):
        for E in (.21,.5,1.,10.):
            s=E*E
            if s<=4*m*m:continue
            forward=annihilation_sigma(s,m)
            reverse=inverse_ee_sigma(s,m)
            assert math.isclose(reverse/forward,
                                1-4*m*m/s,rel_tol=1e-12)
            # Lorentz-invariant thermal 2->2 reaction density has
            # σ(s)*(s-4m²)√s; reverse massless e has σinverse*s√s.
            assert math.isclose(reverse*s,
                                forward*(s-4*m*m),rel_tol=2e-12)


def test_resonance_theta_integral_matches_independent_MB_derivative_for_heat():
    for T in (.5,1.,2.):
        delta=3e-4*T
        plus=s_channel_thermodynamic_MB(T+delta,.02)
        minus=s_channel_thermodynamic_MB(T-delta,.02)
        center=s_channel_thermodynamic_MB(T,.02)
        # Q_forward=−dγ/dβ =T² dγ/dT, independently differentiate
        # the FULL MB reaction density and ensure resonance handling.
        from_fd=T*T*(plus["gamma_forward_MB_GeV4"]-
                      minus["gamma_forward_MB_GeV4"])/(2*delta)
        assert math.isclose(center["energy_forward_MB_GeV5"],
                            from_fd,rel_tol=4e-5),(
            T,center["energy_forward_MB_GeV5"],from_fd)
        assert center["mean_lab_total_energy_per_event_GeV"]>0


def test_forward_dark_vector_resonance_vs_naive_contact_MB_has_huge_real_difference():
    m=.02;T=1.
    r=s_channel_thermodynamic_MB(T,m)
    gprod2=(.1*1e-5*math.sqrt(4*math.pi*ALPHA_EM))**2
    # Independent direct dimensionless-u integral with NO pole,
    #σ_contact = gprod²(s+2m²)/(12π mV⁴ β).
    x=m/T;mV=10.
    def integrand(u):
        s=(T*u)**2
        if s<=4*m*m:return 0.
        beta=math.sqrt(1-4*m*m/s)
        sigma0=gprod2/(12*math.pi)*(s+2*m*m)/(beta*mV**4)
        return sigma0*u*u*(u*u-4*x*x)*kve(1,u)*math.exp(-u)
    cont=quad(integrand,2*x,70,epsabs=0,epsrel=3e-10)[0]
    cont*=4*T**6/(16*math.pi**4)
    ratio=r["gamma_forward_MB_GeV4"]/cont
    print("VECTOR_RESONANT_MB_NEGATIVE_CONTACT",repr(dict(
        gamma_forward=r["gamma_forward_MB_GeV4"],
        gamma_contact=cont,enhancement=ratio,
        forward_over_H_NOT_net=r["forward_over_H_NOT_NET"],
        width_minimal=r["GammaV_minimal_channels_GeV"],
        avg_lab_pair_energy=r["mean_lab_total_energy_per_event_GeV"])),flush=True)
    assert 1400<ratio<1500
    assert math.isclose(r["forward_over_H_NOT_NET"],
                        1593.812,rel_tol=4e-5)
    assert math.isclose(r["gamma_forward_MB_GeV4"],
                        3.450408e-16,rel_tol=5e-5)


def test_BreitWigner_width_controls_resonant_MB_source_as_a_real_unfixed_input():
    w=minimal_vector_width()["width_GeV"]
    a=s_channel_thermodynamic_MB(1.,.02,GammaV=w/2)["gamma_forward_MB_GeV4"]
    b=s_channel_thermodynamic_MB(1.,.02,GammaV=w)["gamma_forward_MB_GeV4"]
    c=s_channel_thermodynamic_MB(1.,.02,GammaV=w*2)["gamma_forward_MB_GeV4"]
    assert a>b>c
    assert 1.6<a/b<2.3
    assert 1.6<b/c<2.3


def test_MB_net_two_temperature_heat_and_entropic_chemical_detailed_balance():
    eq=net_MB_two_temperature_ledger(1.,1.,.02)
    assert eq["R_forward_minus_reverse_GeV4"]==0
    assert eq["Q_dark_energy_loss_GeV5"]==0
    assert eq["H_ideal_MB_GeV4"]==0
    neq=net_MB_two_temperature_ledger(.8,1.,.02)
    assert neq["R_forward_minus_reverse_GeV4"]<0
    assert neq["Q_dark_energy_loss_GeV5"]<0
    assert neq["Q_dark_GeV5"]>0
    assert neq["Q_e_GeV5"]<0
    assert neq["H_ideal_MB_GeV4"]>0
    assert neq["U1_net_defect"]==0
    assert neq["electric_charge_net_defect"]==0
    # Nonzero chemical affinities (µ_chi+µ_antichi) do not assume a
    # single positive-count dark charge across SM pair creation.
    for mu1,mu2 in ((.03,.03),(.03,-.03),(-.1,.02)):
        sameT=net_MB_two_temperature_ledger(
            1.,1.,.02,mu_chi_sum=mu1,mu_e_sum=mu2)
        assert sameT["H_ideal_MB_GeV4"]>=-1e-24
        if mu1==mu2:
            assert abs(sameT["R_forward_minus_reverse_GeV4"])<1e-27
    print("VECTOR_RESONANT_MB_NET_LEDGER",
          repr(dict(R=neq["R_forward_minus_reverse_GeV4"],
                    Qdark=neq["Q_dark_GeV5"],H=neq["H_ideal_MB_GeV4"])),
          flush=True)


def test_resonant_forward_MB_per_chi_T_and_H_gate_not_full_equilibration():
    for T,g,low,hi in ((.1,17.25,.002,.003),
                        (.5,61.75,.9,1.1),
                        (1.,61.75,1500,1700),
                        (2.,75,18000,21000)):
        z=s_channel_thermodynamic_MB(T,.02,gstar=g)
        assert low<z["forward_over_H_NOT_NET"]<hi
        assert z["gamma_forward_MB_GeV4"]>0
        assert z["energy_forward_MB_GeV5"]>0
        assert z["num_theta_segments"]>=10
