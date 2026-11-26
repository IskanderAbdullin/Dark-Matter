"""On-shell Pauli event ledger: independent equilibrium, positivity and recoil."""
import math
import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from microphysics_elastic_event import (
    rest_chi_event_kinematics,elastic_fermion_event)


@pytest.mark.parametrize("z", [0.05,.2,.6,1.])
def test_kinematics_four_momentum_recoil_and_shell(z):
    d=rest_chi_event_kinematics(1.,.2,z)
    assert abs(d["energy_sum_residual"])<2e-15
    assert abs(d["mother_mass_shell_residual"])<3e-15
    assert -1<=d["electron_cos_theta_lab"]<=1
    assert d["chi_recoil_energy"]>0


@pytest.mark.parametrize("T", [.12,.3,1.])
def test_common_temperature_chemical_equilibrium_is_exact(T):
    d=rest_chi_event_kinematics(1.,.2,.5)
    r=elastic_fermion_event(d["chi_before_E"],d["electron_before_E"],
                            d["chi_after_E"],d["electron_after_E"],
                            T,T,muchi=.1,mue=-.04)
    assert abs(r["net_forward_factor"])<2e-15
    assert abs(r["entropy_production_unweighted"])<2e-15


def test_colder_chi_is_heated_by_hotter_e_and_entropy_is_nonnegative():
    d=rest_chi_event_kinematics(1.,.3,.3)
    r=elastic_fermion_event(d["chi_before_E"],d["electron_before_E"],
                            d["chi_after_E"],d["electron_after_E"],
                            .1,.4,muchi=.02,mue=.01)
    assert r["net_forward_factor"]>0
    assert r["chi_energy_transfer_unweighted"]>0
    assert r["electron_energy_transfer_unweighted"]<0
    assert r["entropy_production_unweighted"]>0
    assert abs(r["chi_energy_transfer_unweighted"]+
               r["electron_energy_transfer_unweighted"])<1e-20
    assert math.isclose(
        r["forward_factor"]/r["reverse_factor"],
        math.exp(r["expected_log_forward_reverse"]),rel_tol=1e-12)


def test_hotter_chi_reverses_direction_and_entropy_stays_positive():
    d=rest_chi_event_kinematics(1.,.3,.3)
    r=elastic_fermion_event(d["chi_before_E"],d["electron_before_E"],
                            d["chi_after_E"],d["electron_after_E"],.4,.1)
    assert r["net_forward_factor"]<0
    assert r["chi_energy_transfer_unweighted"]<0
    assert r["entropy_production_unweighted"]>0


def test_reverse_orientation_changes_sign_of_event_source():
    d=rest_chi_event_kinematics(1.,.2,.4)
    a=elastic_fermion_event(d["chi_before_E"],d["electron_before_E"],
                            d["chi_after_E"],d["electron_after_E"],.15,.3)
    rev=elastic_fermion_event(d["chi_after_E"],d["electron_after_E"],
                              d["chi_before_E"],d["electron_before_E"],.15,.3)
    assert math.isclose(a["net_forward_factor"],-rev["net_forward_factor"],
                        rel_tol=1e-12,abs_tol=1e-18)
    assert math.isclose(a["entropy_production_unweighted"],
                        rev["entropy_production_unweighted"],rel_tol=1e-12)


def test_off_shell_and_negative_inputs_rejected():
    with pytest.raises(ValueError):
        elastic_fermion_event(1.,.2,1.01,.1,.2,.3)
    with pytest.raises(ValueError):
        elastic_fermion_event(1.,.2,1.1,.1,0.,.3)
    with pytest.raises(ValueError):
        rest_chi_event_kinematics(1.,.2,1.1)
