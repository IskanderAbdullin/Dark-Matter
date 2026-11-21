"""Finite charge-neutral e± thermal Schur capacity; fixed-nQ not fixed-eta."""
import math
import pytest
from kinetics.fixed_charge_pair_capacity import (
    charge_constrained_heat, independent_finite_difference)
from kinetics.em_pair_thermodynamics import M_E_GEV, equilibrium

@pytest.mark.parametrize("r",[.05,.1,1.])
def test_positive_charge_constrained_heat_and_finite_difference(r):
    T=r*M_E_GEV
    c=charge_constrained_heat(T)
    fd=independent_finite_difference(T)
    assert c["electron_pair_C_at_fixed_net_charge_over_T3"]>0
    assert 0<c["electron_capacity_ratio_CQ_to_Cy"]<=1+1e-12
    assert math.isclose(c["electron_pair_C_at_fixed_net_charge_over_T3"],
                        fd,rel_tol=3e-6)

def test_low_T_charge_constraint_changes_electron_heat_by_orders_of_magnitude():
    r=.025
    c=charge_constrained_heat(r*M_E_GEV)
    fd=independent_finite_difference(r*M_E_GEV)
    assert c["mu_e_minus_over_T"]>10
    assert c["electron_capacity_ratio_CQ_to_Cy"]<.002
    assert math.isclose(c["electron_pair_C_at_fixed_net_charge_over_T3"],
                        fd,rel_tol=3e-6)

def test_near_zero_charge_heat_matches_mu_zero_reference():
    T=M_E_GEV
    c=charge_constrained_heat(T,eta=6e-10)
    assert math.isclose(c["electron_pair_C_at_fixed_net_charge_over_T3"],
                        equilibrium(T)["C_pair_at_fixed_mu0_over_T3"],
                        rel_tol=1e-12)
