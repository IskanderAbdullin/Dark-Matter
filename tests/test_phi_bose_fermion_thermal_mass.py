"""Independent determinants and one-loop T/m threshold for real 86eV scalar.

No claim of a real pole mass or a cosmological symmetry-breaking phase
transition without scalar self-interaction, portals and kinetic history.
"""
import math
import numpy as np
from kinetics.phi_bose_fermion_thermal_mass import (
    M1,M2,MPHI,Y,fermi_I,fermion_mass_eigenvalues,
    fermion_static_curvature,fermion_thermal_potential,total_static_curvature,
    quartic_needed_for_nonnegative_static_curvature,
    critical_curvature_roots,critical_quartic_envelope,
    naive_phi_only_highT_decay_closure,
)


def test_offdiagonal_two_Dirac_mass_matrix_has_exact_quadratic_trace_and_Z2():
    for phi in (.01,.1,.4,1.):
        ev=fermion_mass_eigenvalues(phi)
        rev=fermion_mass_eigenvalues(-phi)
        assert np.allclose(ev,rev,rtol=0,atol=1e-15)
        assert math.isclose(sum(ev**2),M1*M1+M2*M2+2*Y*Y*phi*phi,
                            rel_tol=2e-15)
    assert math.isclose(fermi_I(0.,1.),1/24,rel_tol=5e-13)


def test_analytic_thermal_static_mass_matches_independent_fermion_determinant_second_difference():
    for T in (.01,.02,.035,.04,.1,1.):
        phi=.2
        Vp=fermion_thermal_potential(+phi,T)
        Vm=fermion_thermal_potential(-phi,T)
        V0=fermion_thermal_potential(0.,T)
        assert math.isclose(Vp,Vm,rel_tol=1e-13)
        derivative=(Vp+Vm-2*V0)/phi**2
        exact=fermion_static_curvature(T)
        assert exact!=0
        assert math.isclose(derivative,exact,rel_tol=8e-4)


def test_finite_mass_negative_static_curvature_and_tiny_vacuum_scalar_mass_gate():
    assert total_static_curvature(.002)>0
    for T in (.005,.01,.02,.03,.035):
        assert total_static_curvature(T)<0
    for T in (.04,.1,1.):
        assert total_static_curvature(T)>0
    low,high=critical_curvature_roots()
    assert math.isclose(low,.0032453669477629813,rel_tol=2e-8)
    assert math.isclose(high,.03548213975072978,rel_tol=2e-8)
    assert math.isclose(fermion_static_curvature(.02),
                        -1.7913210999114e-11,rel_tol=2e-8)


def test_highT_two_Dirac_coefficient_and_static_scale_much_greater_than_86eV():
    for T,relative in ((1.,.02),(10.,3e-4),(100.,2e-5)):
        d=fermion_static_curvature(T)
        assert abs(d/(Y*Y*T*T/3)-1)<relative
    eff=math.sqrt(total_static_curvature(1.))
    assert math.isclose(eff,0.000573685,rel_tol=3e-6)
    assert eff/MPHI>6000
    # The following is a phi-only hypothetical substitution, NOT an
    # actual full thermal kinematic closing of the 1->2 reaction:
    assert math.isclose(naive_phi_only_highT_decay_closure(),
                        138.5640646054301,rel_tol=1e-12)


def test_scalar_quartic_missing_parameter_controls_conditional_negative_window():
    env=critical_quartic_envelope()
    assert .017<env["temperature_GeV"]<.020
    assert math.isclose(env["min_quartic_lambda_phi4_over_4factorial"],
                        1.0851571368922e-6,rel_tol=2e-7)
    for T in np.geomspace(.002,1.,40):
        assert total_static_curvature(T,lam=1.1e-6)>0
    assert quartic_needed_for_nonnegative_static_curvature(.02)>1e-6
