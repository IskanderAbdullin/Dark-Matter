"""Independent gamma-trace / thermal MB checks of off-diagonal Yukawa 2->2.

A valid vacuum pair-annihilation channel is not proof of net thermal
chemical equilibration: backreaction/Bose/Pauli and SM/FRW absent.
"""
import math
import numpy as np
from scipy.integrate import quad
from scipy.special import kv
from kinetics.phi_bose_pair_annihilation_mb import (
    M1,M2,MPHI,Y,dot,scattering_vectors,spin_averaged_amp2,
    cross_section,thermal_forward_rate_MB,
)


def _gamma_matrices():
    z=np.zeros((2,2),complex)
    i2=np.eye(2,dtype=complex)
    pauli=[np.array([[0,1],[1,0]],complex),
           np.array([[0,-1j],[1j,0]],complex),
           np.diag([1,-1]).astype(complex)]
    gam0=np.block([[i2,z],[z,-i2]])
    gs=[gam0]+[np.block([[z,s],[-s,z]]) for s in pauli]
    return gs


def _gamma_trace_independent(s,c):
    """Direct 4x4 Dirac trace, not relying on covariant dot trace formula."""
    p1,p2,k1,k2=scattering_vectors(s,c)
    gs=_gamma_matrices()
    def slash(v):
        return gs[0]*v[0]-sum((gs[j]*v[j] for j in (1,2,3)),
                              start=np.zeros((4,4),complex))
    I=np.eye(4,dtype=complex)
    q1=p1-k1
    q2=p1-k2
    A=(slash(q1)+M1*I)/(dot(q1,q1)-M1*M1)+(
       slash(q2)+M1*I)/(dot(q2,q2)-M1*M1)
    # Dirac adjoint of A is also A since all real scalar denominators:
    # Abar=gamma0 A† gamma0; verify separately as a source sign gate.
    Abar=gs[0]@A.conj().T@gs[0]
    return float((Y**4/4*np.trace(
        (slash(p2)-M2*I)@A@(slash(p1)+M2*I)@Abar)).real)


def test_invariant_trace_equals_independent_gamma_matrix_spin_sum():
    for E in (.041,.05,.1,.4,1.,2.):
        for c in (0.,.3,.9):
            s=E*E
            got=spin_averaged_amp2(s,c)
            reference=_gamma_trace_independent(s,c)
            assert got>0 and reference>0
            assert math.isclose(got,reference,rel_tol=2e-12)


def test_identical_real_scalar_half_factor_and_angle_quadrature():
    for E,expected in ((.05,1.4462000244688534e-12),
                       (.1,1.2720976042160392e-12),
                       (.5,8.584263024698827e-14),
                       (1.,2.935350794935496e-14)):
        coarse=cross_section(E*E,nangle=96)
        fine=cross_section(E*E,nangle=128)
        assert fine>0
        assert math.isclose(fine,expected,rel_tol=1e-8)
        assert math.isclose(coarse,fine,rel_tol=2e-7)
    assert math.isclose(cross_section(1.,Y/10)/cross_section(1.,Y),
                        1e-4,rel_tol=2e-12)
    assert cross_section((2*M2)**2)==0


def test_exact_Gondolo_Gelmini_MB_reaction_density_vs_separate_sigma_v_integral():
    T=1.
    x=M2/T
    # Independent thermal reaction formula in center of mass E=√s,
    # not reusing the production module's dimensionless-u integral.
    integral=quad(
       lambda E: 2*E*cross_section(E*E)*(E*E-4*M2*M2)*E*kv(1,E/T),
       2*M2,40*T,epsabs=1e-35,epsrel=5e-6,limit=140)[0]
    gamma=(2*2)*T/(32*math.pi**4)*integral
    ref=thermal_forward_rate_MB(1.,61.75)
    assert math.isclose(gamma,ref["forward_pair_reaction_density_GeV4"],
                        rel_tol=5e-5)
    assert math.isclose(ref["forward_pair_rate_per_chi2_GeV"],
                        1.0308030420839305e-15,rel_tol=4e-5)
    assert math.isclose(ref["forward_pair_rate_over_H"],
                        964.7820534903847,rel_tol=4e-5)


def test_MB_pair_forward_rate_exceeds_H_for_strong_y_but_not_weak_y():
    for T,gstar,lo,hi in ((.1,17.25,6900,7000),
                           (1.,61.75,950,980),
                           (100.,106.75,23,24)):
        r=thermal_forward_rate_MB(T,gstar,nangle=96)
        assert lo<r["forward_pair_rate_over_H"]<hi
        assert r["forward_pair_reaction_density_GeV4"]>0
    # Scales y^4 in forward vacuum cross-section and thermal rate.
    strong=thermal_forward_rate_MB(1.,61.75,y=Y)
    feeble=thermal_forward_rate_MB(1.,61.75,y=Y*1e-6)
    assert math.isclose(feeble["forward_pair_rate_over_H"]/
                        strong["forward_pair_rate_over_H"],1e-24,
                        rel_tol=2e-9)
    assert feeble["forward_pair_rate_over_H"]<1e-18


def test_raw_cosine_gauss_underresolves_ultrarelativistic_t_channel_peaks():
    # Deliberately reproduce and reject the old 36-node plain-cos(theta)
    # estimate: it undercounts 100-GeV scattering by ~46%, unlike
    # logarithmically clustered hyperbolic-angle integration.
    from numpy.polynomial.legendre import leggauss
    z,w=leggauss(36)
    for E,max_fraction in ((5.,.96),(100.,.60)):
        s=E*E
        pin=math.sqrt(s/4-M2*M2)
        pout=math.sqrt(s/4-MPHI*MPHI)
        old=sum(wi*spin_averaged_amp2(s,float(c)) for c,wi in zip(z,w))
        old*=pout/(64*math.pi*s*pin)
        new=cross_section(s,nangle=128)
        assert old/new<max_fraction
    assert math.isclose(cross_section(100.**2,nangle=96),
                        cross_section(100.**2,nangle=128),rel_tol=1e-5)
