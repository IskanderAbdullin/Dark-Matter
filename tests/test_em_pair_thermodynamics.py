"""Independent massive e±-photon EOS and neutrality benchmarks."""
import math
import pytest
from kinetics.em_pair_thermodynamics import (
    M_E_GEV, equilibrium, bessel_pair, charge_neutral_mu, N_PHOTON_OVER_T3)

@pytest.mark.parametrize("ratio",[.1,1.,10.])
def test_independent_bessel_series_and_capacity_derivative(ratio):
    m=M_E_GEV
    row=equilibrium(ratio*m)
    bp,br=bessel_pair(1/ratio)
    assert math.isclose(row["P_pair_over_T4"],bp,rel_tol=3e-11)
    assert math.isclose(row["rho_pair_over_T4"],br,rel_tol=3e-11)
    h=1e-4
    Tp=ratio*m*(1+h);Tm=ratio*m*(1-h)
    rhop=equilibrium(Tp)["rho_pair_over_T4"]*Tp**4
    rhom=equilibrium(Tm)["rho_pair_over_T4"]*Tm**4
    fd=(rhop-rhom)/(Tp-Tm)/(ratio*m)**3
    assert math.isclose(fd,row["C_pair_at_fixed_mu0_over_T3"],rel_tol=5e-7)

def test_high_and_low_temperature_gstar_limits():
    assert abs(equilibrium(100*M_E_GEV)["gstar_rho_gamma_plus_pair"]-5.5)<5e-5
    assert abs(equilibrium(.05*M_E_GEV)["gstar_rho_gamma_plus_pair"]-2)<5e-5
    for r in [.05,.1,1.,10.]:
        a=equilibrium(r*M_E_GEV)
        assert a["C_gamma_plus_pair_over_T3"]>0
        assert a["gstar_s_gamma_plus_pair"]>0

def test_charge_neutrality_requires_nonzero_electron_mu_below_annihilation():
    eta=6e-10
    warm=charge_neutral_mu(.1*M_E_GEV,eta=eta)
    cold=charge_neutral_mu(.025*M_E_GEV,eta=eta)
    for a in (warm,cold):
        assert math.isclose(a["net_electron_charge_over_T3"],
                            eta*N_PHOTON_OVER_T3,rel_tol=1e-9)
    assert warm["mu_e_minus_over_T"]<1e-5
    assert cold["mu_e_minus_over_T"]>10.
    assert cold["positron_to_electron_density_ratio"]<1e-10
