"""Independent vacuum t/u gamma-matrix trace and conditional stable χ2 relic.

Strictly tests *restricted* symmetric thermal MB + fixed-g* p-wave
approximation, not the full dark+SM FRW, BEC/Ω_phi/CLASS-ready.
"""
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from kinetics.phi_bose_stable_chi2_relic_gate import (
    exact_scalar_pwave_b,vector_ee_swave_a,conditional_symmetric_relic,
    coupling_for_thermal_relic_target,mb_yeq_single_dirac,
)
from kinetics.phi_bose_pair_independent_vectorized import sigma_independent
from kinetics.phi_bose_vector_portal_annihilation_resonance import (
    annihilation_sigma,
)


def _gamma():
    I=np.eye(2,dtype=complex);Z=np.zeros((2,2),complex)
    sig=(np.array([[0,1],[1,0]],complex),
         np.array([[0,-1j],[1j,0]],complex),
         np.diag([1,-1]).astype(complex))
    return [np.block([[I,Z],[Z,-I]])]+[
        np.block([[Z,x],[-x,Z]]) for x in sig]


def _direct_gamma_trace_sigma_v(beta,m=.02,M=.1,mphi=86e-9,y=.001):
    """Independently spin-average the FULL four-dimensional Dirac trace."""
    gam=_gamma();I=np.eye(4,dtype=complex)
    def slash(q):
        return sum((factor*gam[i]*q[i] for i,factor in
                   enumerate((1,-1,-1,-1))),start=np.zeros((4,4),complex))
    def mink(x,z):
        return float(x[0]*z[0]-np.dot(x[1:],z[1:]))
    E=m/math.sqrt(1-beta*beta)
    p=beta*E
    k=math.sqrt(E*E-mphi*mphi)
    parent=np.array([E,0,0,p])
    anti=np.array([E,0,0,-p])
    x,w=leggauss(64)
    integral=0.
    for c,wc in zip(x,w):
        trans=k*math.sqrt(1-c*c)
        products=(np.array([E,trans,0,k*c]),
                  np.array([E,-trans,0,-k*c]))
        exchange=np.zeros((4,4),complex)
        for daughter in products:
            q=parent-daughter
            exchange+=(slash(q)+M*I)/(mink(q,q)-M*M)
        amp=y**4*float(np.trace((slash(anti)-m*I)@exchange@
                  (slash(parent)+m*I)@exchange).real/4)
        assert amp>=-1e-20
        integral+=wc*amp
    # Identical final real scalar 1/2 included.
    return integral*k/(32*math.pi*4*E*E*E)


def test_analytic_pwave_from_independent_complex_gamma_4x4_trace():
    exact=exact_scalar_pwave_b()
    assert math.isclose(exact,1.3130827280563007e-12,rel_tol=3e-13)
    for beta in (.2,.1,.05,.02):
        vM=2*beta
        direct=_direct_gamma_trace_sigma_v(beta)
        old=2*beta*sigma_independent(
            4*(.02/math.sqrt(1-beta*beta))**2,angle_order=128)
        assert math.isclose(direct,old,rel_tol=2e-10),(beta,direct,old)
        assert abs(direct/(exact*vM*vM)-1)<.18*beta*beta+1e-8
    # Extrapolate p->0 as physical p-wave, not constant s-wave.
    assert math.isclose(
        _direct_gamma_trace_sigma_v(.02)/(.04**2),exact,rel_tol=1e-4)


def test_vector_e_cross_section_true_swave_is_negligible_to_scalar_pwave_at_x8():
    a=vector_ee_swave_a()
    b=exact_scalar_pwave_b()
    assert math.isclose(a,1.1676134457946275e-21,rel_tol=5e-11)
    for beta in (.1,.03,.01):
        m=.02
        s=4*m*m/(1-beta*beta)
        sig=2*beta*annihilation_sigma(s,m)
        assert math.isclose(sig,a,rel_tol=.6*beta*beta+1e-6)
    assert a/(6*b/8)<1e-8


def test_thermal_MB_per_charge_density_and_relic_ode_initial_conditions():
    y3=mb_yeq_single_dirac(3.)
    assert 0<y3<.02
    assert mb_yeq_single_dirac(6.)<y3
    a=conditional_symmetric_relic(x_initial=3.)
    b=conditional_symmetric_relic(x_initial=5.)
    assert math.isclose(a["Y_chi2_per_charge_final"],
                        b["Y_chi2_per_charge_final"],rel_tol=3e-8)
    assert a["Y_chi2_per_charge_final"]>0
    assert math.isclose(a["Omega_chi2_plus_anti2_h2"],761.365,
                        rel_tol=5e-5)
    assert a["relative_to_OmegaDM_0p120"]>6000
    assert not a["class_ready"] and not a["physical_thermal_history_solved"]


def test_late_freezeout_yield_converges_and_symmetric_chi2_anti2_counted_twice():
    a=conditional_symmetric_relic(x_final=10000.)
    b=conditional_symmetric_relic(x_final=100000.)
    assert abs(a["Omega_chi2_plus_anti2_h2"]/b[
        "Omega_chi2_plus_anti2_h2"]-1)<1e-6
    factor=2891.2/1.05375e-5
    assert math.isclose(b["Omega_chi2_plus_anti2_h2"],
                        2*.02*factor*b["Y_chi2_per_charge_final"],
                        rel_tol=1e-14)


def test_chi2_coupling_scan_boltzmann_relic_decreases_without_changing_86eV_mass():
    r=[conditional_symmetric_relic(y=y) for y in (.001,.004,.01,.014,.02)]
    omega=[x["Omega_chi2_plus_anti2_h2"] for x in r]
    assert all(a>b for a,b in zip(omega,omega[1:]))
    assert math.isclose(omega[0],761.365,rel_tol=5e-5)
    assert .09<omega[3]<.11
    # A different coupling could evade THIS toy overclosure,
    # but changes all Yukawa decay+pair rates and phi production.
    target=coupling_for_thermal_relic_target()
    assert .0132<target<.0134
    assert math.isclose(conditional_symmetric_relic(y=target)[
        "Omega_chi2_plus_anti2_h2"],.12,rel_tol=1e-4)


def test_conditional_cosmological_assumptions_are_explicit():
    r=conditional_symmetric_relic()
    assert r["gstar_rho_assumed"]==10.75
    assert r["gstar_S_assumed"]==10.75
    assert r["mphi_GeV"]==86e-9
    assert not r["physical_thermal_history_solved"]
    # This test deliberately does NOT infer a full thermal chi2 population
    # from a forward/H event rate at GeV, and does not assert any
    # actual BEC, phi Omega, late frozen f(q) or observable fit.
