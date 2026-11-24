"""Independent test of chosen vacuum Yukawa y against RD Hubble / old FRW κ.

Neither the benchmark width nor MB test-particle rate establishes physical
full χ-SM in-medium reaction, a successful freeze-in model or dark relic.
"""
import math
from scipy.integrate import quad
from kinetics.phi_bose_physical_cutcell_weak import WIDTH
from kinetics.phi_bose_yukawa_hubble_gate import (
    M1,M2,MPHI,YBENCH,rest_width,hubble_rd,rate_gate,
    mb_mean_rest_time_dilation,mb_crossover_T,
)


def test_spin_averaged_dirac_scalar_tree_width_matches_independent_weak_source():
    gamma=rest_width()
    assert math.isclose(gamma,WIDTH,rel_tol=1e-14)
    assert math.isclose(gamma,2.7501974166242436e-9,rel_tol=1e-13)
    assert math.isclose(rest_width(1e-4)/gamma,.01,rel_tol=1e-14)
    assert math.isclose(rest_width(0.),0.,abs_tol=0.)
    # Independent rest-frame spin trace:
    # Spin-averaged |M|² = y²[(M1+M2)²−mphi²].
    amplitude_sq=YBENCH**2*((M1+M2)**2-MPHI**2)
    kstar=math.sqrt(
        (M1*M1-(M2+MPHI)**2)*(M1*M1-(M2-MPHI)**2)
    )/(2*M1)
    assert math.isclose(gamma,amplitude_sq*kstar/(8*math.pi*M1*M1),
                        rel_tol=2e-14)


def test_exact_MBoltzmann_K1K2_time_dilation_vs_direct_integrals():
    # q=p/T, x=M1/T: compare independently evaluated physical p²dp
    # MB number and time-dilated number integrals with analytic K1/K2.
    for T in (.1,1.,100.):
        x=M1/T
        n=quad(lambda q:q*q*math.exp(-math.hypot(q,x)),
               0.,80.,epsabs=1e-13,epsrel=1e-12)[0]
        weighted=quad(lambda q:q*q*x/math.hypot(q,x)*
                       math.exp(-math.hypot(q,x)),
                       0.,80.,epsabs=1e-13,epsrel=1e-12)[0]
        assert math.isclose(weighted/n,
                            mb_mean_rest_time_dilation(T),
                            rel_tol=2e-10)


def test_ybenchmark_is_NOT_prior_small_kappa_FRW_pilot():
    for gstar in (17.25,61.75):
        row=rate_gate(.1,gstar)
        assert row["Gamma_rest_over_H"]>2e11
        assert row["Gamma_MB_per_parent_over_H"]>9e10
        assert row["benchmark_y_for_kappa_rest_0p1"]<1e-9
        assert math.isclose(
            rate_gate(.1,gstar,row["benchmark_y_for_kappa_rest_0p1"])
                ["Gamma_rest_over_H"],.1,rel_tol=1e-13)
    g1=rate_gate(1.,61.75)
    assert math.isclose(g1["Gamma_rest_over_H"],
                        2574052464.7177434,rel_tol=1e-12)
    assert math.isclose(g1["Gamma_MB_per_parent_over_H"],
                        127136889.18737303,rel_tol=1e-12)


def test_high_temperature_MB_test_parent_rate_and_crossover_gstar_fixed():
    t=mb_crossover_T(106.75)
    assert math.isclose(t,460.8652365658074,rel_tol=1e-10)
    ratio=rate_gate(t,106.75)["Gamma_MB_per_parent_over_H"]
    assert math.isclose(ratio,1.,rel_tol=1e-11)
    assert rate_gate(100.,106.75)["Gamma_MB_per_parent_over_H"]>90
    assert rate_gate(1000.,106.75)["Gamma_MB_per_parent_over_H"]<.1
    # MB time dilation is parametrically M/(2T) at large T; do not
    # equate Γ_rest/H with Γ_lab/H in a relativistic parent plasma.
    for T in (100.,500.,1000.):
        assert math.isclose(mb_mean_rest_time_dilation(T),M1/(2*T),
                            rel_tol=5e-6)
