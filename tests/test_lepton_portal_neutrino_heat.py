"""Chiral neutrino vector portal and independent 5D thermal heat controls."""
import math
import pytest
from kinetics.lepton_portal_neutrino_heat import (
    massive_chi_neutrino_trace,heat_lepton_comparison)

@pytest.mark.parametrize("c",[-.8,-.2,.3,.9])
def test_single_left_helicity_trace_equals_dirac_vector_average_with_vector_chi(c):
    V,L,analytic,imag=massive_chi_neutrino_trace(.1,.5,c)
    assert imag<1e-12
    assert math.isclose(V,L,rel_tol=2e-13)
    assert math.isclose(V,analytic,rel_tol=2e-13)

@pytest.mark.parametrize("ratio,lower,upper",[(.1,1000,1200),(1,1.1,1.3),(10,1,1.01)])
def test_conditional_two_nu_flavor_vs_e_heat_without_dark_relic_claim(ratio,lower,upper):
    a=heat_lepton_comparison(ratio,n=20)
    assert a["Q_chi_two_active_nu_G2_me9_scaled"]>0
    assert a["Q_chi_e_minus_e_plus_G2_me9_scaled"]>0
    assert lower<a["nu_over_e_heat_ratio"]<upper

def test_heat_signs_near_T_eq_in_both_lepton_channels():
    a=heat_lepton_comparison(.1,Tchi_over_T=1.1,n=20)
    assert a["Q_chi_two_active_nu_G2_me9_scaled"]<0
    assert a["Q_chi_e_minus_e_plus_G2_me9_scaled"]<0
