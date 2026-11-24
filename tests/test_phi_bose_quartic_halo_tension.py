"""Independent identical-scalar tree cross-section and conditional halo bound."""
import math
from kinetics.phi_bose_fermion_thermal_mass import (
    MPHI,critical_quartic_envelope,total_static_curvature
)
from kinetics.phi_bose_quartic_halo_tension import (
    CM2_G_PER_GEV_MINUS3,sigma_contact_GeVminus2,sigma_over_m_cm2_g,
    lambda_for_reference_sigma_m,fixed_lambda_mass_for_reference_sigma_m,audit
)

def test_identical_real_scalar_half_factor_and_units():
    lam=1e-6
    s=4*MPHI**2
    independent=4*math.pi*lam**2/(128*math.pi**2*s)
    assert math.isclose(independent,sigma_contact_GeVminus2(lam),rel_tol=1e-14)
    assert math.isclose(CM2_G_PER_GEV_MINUS3,2.1844e-4,rel_tol=5e-4)
    assert math.isclose(sigma_over_m_cm2_g(lam),853.9820784064723,rel_tol=3e-12)

def test_fermion_thermal_stabilizing_quartic_has_large_Born_cluster_sigma():
    lam=critical_quartic_envelope()["min_quartic_lambda_phi4_over_4factorial"]
    assert math.isclose(lam,1.0851571368922225e-6,rel_tol=3e-7)
    assert math.isclose(sigma_over_m_cm2_g(lam),
                        1005.620270173485,rel_tol=7e-7)
    assert total_static_curvature(.02,lam=lam)>0

def test_cluster_reference_allowed_lambda_is_insufficient_for_conditional_20MeV_Z2_static_gate():
    lam=critical_quartic_envelope()["min_quartic_lambda_phi4_over_4factorial"]
    for upper,expected in ((1.,3.421965460115207e-8),
                           (2.,4.839389963667214e-8),
                           (5.,7.651747385473947e-8)):
        maxlam=lambda_for_reference_sigma_m(upper)
        assert math.isclose(maxlam,expected,rel_tol=3e-12)
        assert math.isclose(sigma_over_m_cm2_g(maxlam),upper,rel_tol=2e-12)
        assert lam/maxlam>14
        assert total_static_curvature(.02,lam=maxlam)<0

def test_scalar_mass_cubed_and_quartic_squared_scaling():
    lam=1e-6
    for scale in (2.,4.,8.):
        assert math.isclose(sigma_over_m_cm2_g(lam,MPHI*scale),
                            sigma_over_m_cm2_g(lam)/scale**3,rel_tol=2e-14)
    x=fixed_lambda_mass_for_reference_sigma_m(
        critical_quartic_envelope()["min_quartic_lambda_phi4_over_4factorial"],1.)
    assert math.isclose(x*1e9,861.6081351505869,rel_tol=3e-7)

def test_conditional_report_has_no_physical_relic_claim():
    a=audit()
    assert a["sigma_over_mass_cm2_g"]>1000
    assert a["mphi_vacuum_eV"]==86.
    assert len(a["reference_scenarios"])==4
    assert a["class_ready"] is False
