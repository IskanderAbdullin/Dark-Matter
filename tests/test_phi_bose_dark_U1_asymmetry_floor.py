"""Exact dark U1 theorem and asymmetry relic floor; independent of MB toy."""
import math
import pytest
from kinetics.phi_bose_dark_U1_asymmetry_floor import (
    exactly_conserved_dark_charge,minimum_charged_fermion_omega_h2,
    maximum_asymmetry_given_fraction,verify_exact_event_charge_conservation,
    MCHI2_GEV,S0_CM3,RHOC_OVER_H2_GEV_CM3,
)


def test_exact_noether_integer_stoichiometry_all_real_yukawa_and_vector_channels():
    d=verify_exact_event_charge_conservation()
    assert len(d)==7
    assert all(x==0 for x in d.values())
    assert exactly_conserved_dark_charge([1,0,0,0,2])==1
    assert exactly_conserved_dark_charge([0,0,1,0,3])==1
    assert exactly_conserved_dark_charge([3,2,10,12,300])==-1


def test_asymmetric_charged_partner_lower_bound_by_triangle_inequality():
    # Integer and fractional physical positive comoving populations.
    for a in ((4.,1.,2.,3.,5.),(.1,4.,2.,.3,0.),(0.,0.,0.,0.,25.)):
        Q=exactly_conserved_dark_charge(a)
        energy=(a[0]+a[1])*.1+(a[2]+a[3])*.02
        assert energy+1e-14>=.02*abs(Q)
    # Convert to today's actual entropy normalization. Never assume
    # that early SM comoving entropy was separately constant.
    for eta in (-3e-8,-1e-10,0,1e-10,3e-8):
        omega=minimum_charged_fermion_omega_h2(eta)
        assert omega>=0
        assert math.isclose(omega,
                            MCHI2_GEV*abs(eta)*S0_CM3/
                            RHOC_OVER_H2_GEV_CM3,rel_tol=2e-15)


def test_maximum_asymmetry_for_1_percent_and_full_dark_abundance():
    one=maximum_asymmetry_given_fraction(1.)
    per=maximum_asymmetry_given_fraction(.01)
    assert math.isclose(one,2.1868082457111234e-8,rel_tol=2e-6)
    assert math.isclose(per,one/100,rel_tol=2e-15)
    assert math.isclose(minimum_charged_fermion_omega_h2(one),.12,
                        rel_tol=1e-14)
    assert math.isclose(minimum_charged_fermion_omega_h2(per),.0012,
                        rel_tol=1e-14)
    assert minimum_charged_fermion_omega_h2(1e-9)>.005


def test_strict_scope_does_not_pretend_measured_dark_asymmetry_or_phi_relic():
    with pytest.raises(ValueError):
        exactly_conserved_dark_charge([1.,-1.,1.,1.,0.])
    with pytest.raises(ValueError):
        maximum_asymmetry_given_fraction(1.5)
    with pytest.raises(ValueError):
        minimum_charged_fermion_omega_h2(float("nan"))
