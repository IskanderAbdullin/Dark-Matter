"""Independent Bessel and MB phase-space checks, not quantum mediator plasma."""
import math
import pytest
from kinetics.vector_mediator_inverse_decay_gate import (
 analytic_time_dilation,boltzmann_time_dilation,conditional_neutrino_decay)

@pytest.mark.parametrize("x",[.03,.1,.3,1,3])
def test_independent_MB_time_dilation_Bessel_K1_K2(x):
    a=analytic_time_dilation(x,1)
    b=boltzmann_time_dilation(x,1)
    assert math.isclose(a,b,rel_tol=1e-9)
    assert 0<a<1

@pytest.mark.parametrize("mV",[.03,.1,.3])
def test_fixed_contact_product_perturbative_chi_minimal_gprime_has_large_MB_GammaH(mV):
    a=conditional_neutrino_decay(mV_me=mV)
    assert math.isclose(a["gprime_required"]*a["gchi_effective_assumed"],
                        a["fixed_coupling_product"],rel_tol=1e-14)
    assert a["Gamma_MB_over_H"]>1000
    assert math.isclose(a["existing_MB_vector_Pauli_decay_lower_over_H_at_mu_nu0"],
                        a["Gamma_MB_over_H"]/4,rel_tol=1e-15)
    b=conditional_neutrino_decay(mV_me=mV,gchi_eff=.1)
    assert b["Gamma_MB_over_H"]>a["Gamma_MB_over_H"]
    assert math.isclose(a["Gamma_MB_over_H"]/b["Gamma_MB_over_H"],
                        (.1/math.sqrt(4*math.pi))**2,rel_tol=1e-12)
