"""Conditional finite-vector electron and neutrino split thermal bath audit."""
import math
from kinetics.finite_mediator_split_lepton_reference import (
    split_heat,TNU_OVER_TGAMMA,
)


def test_conditional_split_electron_neutrino_heat_sign():
    a=split_heat(mV=.1,nr=5,nmu=5,nphi=12)
    assert a["Qe"]>0
    assert a["Qtwo_active_nu"]<0
    assert a["Q_total"]<0
    assert math.isclose(a["Q_total"],a["Qe"]+a["Qtwo_active_nu"],rel_tol=1e-15)
    assert math.isclose(a["Tnu_over_Tgamma"],TNU_OVER_TGAMMA)


def test_equal_bath_heat_sign_differs_from_split_bath():
    equal=split_heat(mV=.1,nr=5,nmu=5,nphi=12,Tnu_over_Tgamma=1.)
    assert equal["Qe"]>0 and equal["Qtwo_active_nu"]>0
    assert equal["Q_total"]>0


def test_coarse_nonconstant_transfer_heat_is_not_a_physical_precision_claim():
    a=split_heat(mV=.1,nr=5,nmu=5,nphi=12)
    b=split_heat(mV=.1,nr=5,nmu=5,nphi=12,
                 polar="transfer",nang=24)
    assert a["Q_total"]<0 and b["Q_total"]<0
    assert abs(b["Q_total"]/a["Q_total"]-1)>.01
