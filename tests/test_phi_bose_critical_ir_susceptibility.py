"""Independent physical p-space Gauss vs IR-core stretched Bose susceptibility.

No BEC proof; the first hard-P pair-source cell is collision-inactive.
"""
import math
import numpy as np
import mpmath as mp
from kinetics.phi_bose_physical_cutcell_weak import MPHI
from kinetics.phi_bose_critical_ir_susceptibility import (
    first_frozen_cell,stretched_physical_moments,
    legendre_physical_moments,leading_asymptotic_G00,
)


def test_frozen_ir_first_cell_is_exactly_below_pair_kinematic_source():
    lo,hi=first_frozen_cell()
    assert math.isclose(lo,86e-9,rel_tol=1e-13)
    assert math.isclose(hi,9.999751861499219e-5,rel_tol=1e-12)
    assert hi>lo


def test_18_gauss_number_looks_converged_while_true_bose_G00_is_underresolved():
    d=1e-16
    correct=stretched_physical_moments(d)
    n18=legendre_physical_moments(d,18)
    n576=legendre_physical_moments(d,576)
    assert math.isclose(correct["susceptibility_G00_GeV3"],
                        .00028488368550289877,rel_tol=2e-10)
    assert .99999<n18["first_number_moment_M0_GeV3"] / \
                    correct["first_number_moment_M0_GeV3"]<1.00001
    assert n18["susceptibility_G00_GeV3"] / \
           correct["susceptibility_G00_GeV3"]<.004
    assert n576["susceptibility_G00_GeV3"] / \
           correct["susceptibility_G00_GeV3"]<.04


def test_true_critical_massive_Bose_susceptibility_grows_like_inverse_sqrt_delta():
    a=stretched_physical_moments(1e-16)["susceptibility_G00_GeV3"]
    b=stretched_physical_moments(1e-18)["susceptibility_G00_GeV3"]
    assert 9.9<b/a<10.1
    assert math.isclose(leading_asymptotic_G00(1e-18)/
                        b,1.,rel_tol=5e-4)
    assert 0<b and 0<a


def test_independent_50_digit_mpmath_momentum_transform():
    mp.mp.dps=50
    m=mp.mpf("8.6e-8")
    d=mp.mpf("1e-16")
    hi=mp.mpf(str(first_frozen_cell()[1]))
    pc=mp.sqrt(2*m*d)
    pmax=mp.sqrt((hi-m)*(hi+m))
    umax=mp.asinh(pmax/pc)
    def f(u):
        p=pc*mp.sinh(u)
        e=mp.sqrt(p*p+m*m)
        k=p*p/(e+m)+d
        b=1/mp.expm1(k)
        hat=(1-p*p/(e+m)/(hi-m))**2
        return (p*p*pc*mp.cosh(u)/(2*mp.pi**2))*b*(1+b)*hat*hat
    decimal=mp.quad(f,[0,1,2,4,8,12,16,umax])
    fast=stretched_physical_moments(float(d))["susceptibility_G00_GeV3"]
    assert abs(float(decimal)/fast-1)<3e-10


def test_not_just_more_gauss_nodes_can_fix_nearcritical_core():
    d=1e-18
    exact=stretched_physical_moments(d)["susceptibility_G00_GeV3"]
    ratios=[legendre_physical_moments(d,n)["susceptibility_G00_GeV3"]/exact
            for n in (18,72,288,576)]
    assert all(0<r<.01 for r in ratios)
    assert all(b>a for a,b in zip(ratios,ratios[1:]))


def test_fixed_positive_geometric_quadrature_resolves_all_subcritical_IR_gaps_without_regridding():
    from kinetics.phi_bose_critical_ir_susceptibility import (
        composite_physical_moments,fixed_geometric_momentum_gauss_grid
    )
    p,w=fixed_geometric_momentum_gauss_grid()
    assert len(p)==len(w)==168
    assert np.min(p)>0 and np.min(w)>0
    # Same nodes are used for all chemical potentials; this is required
    # for a self-consistent discrete physical entropy/H/Jacobian identity.
    ratios=[]
    for delta in (MPHI,1e-10,1e-12,1e-14,1e-16,1e-18):
        true=stretched_physical_moments(delta)
        got=composite_physical_moments(delta)
        ratios.append(got["susceptibility_G00_GeV3"]/
                      true["susceptibility_G00_GeV3"])
        assert abs(got["relative_DOS_capacity_error"])<3e-14
        assert abs(got["first_number_moment_M0_GeV3"]/
                   true["first_number_moment_M0_GeV3"]-1)<3e-10
    assert max(abs(x-1) for x in ratios)<2e-7


def test_coarse_geometric_partition_has_measurable_IR_bias_despite_exact_capacity():
    from kinetics.phi_bose_critical_ir_susceptibility import composite_physical_moments
    ref=stretched_physical_moments(1e-16)["susceptibility_G00_GeV3"]
    rough=composite_physical_moments(1e-16,n_geometric=12,nper=6)
    fine=composite_physical_moments(1e-16,n_geometric=20,nper=8)
    assert abs(rough["relative_DOS_capacity_error"])<3e-14
    assert abs(fine["relative_DOS_capacity_error"])<3e-14
    assert 4e-4 < abs(rough["susceptibility_G00_GeV3"]/ref-1)<7e-4
    assert abs(fine["susceptibility_G00_GeV3"]/ref-1)<2e-7
