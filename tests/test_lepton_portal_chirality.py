"""Lepton doublet gauge-completion weight; conditional massless equal-T limit."""
import math
import sys
from pathlib import Path
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from lepton_portal_chirality import same_temperature_massless_lepton_weights
from relativistic_weak_form_2to2 import weak_energy_moments


def test_vector_electron_requires_equal_L_R_and_neutrino_channel():
    r=same_temperature_massless_lepton_weights(1e-3,1e-3)
    assert r["G_e_axial"]==0
    assert math.isclose(r["G_e_vector"],1e-3)
    assert r["electron_only_weight"]==4e-6
    assert r["leptons_completed_weight"]==6e-6
    assert r["ratio_full_over_electron_only"]==1.5


def test_right_only_has_no_neutrino_but_is_not_electron_vector():
    r=same_temperature_massless_lepton_weights(0,1e-3)
    assert r["neutrino_left_weight"]==0
    assert r["ratio_full_over_electron_only"]==1
    assert r["G_e_axial"]!=0


def test_left_only_doubles_weight_conditional_common_bath():
    r=same_temperature_massless_lepton_weights(1e-3,0)
    assert r["electron_left_weight"]==r["neutrino_left_weight"]
    assert r["ratio_full_over_electron_only"]==2


def test_same_T_additional_channel_in_scattering_weak_heat():
    a=weak_energy_moments(M=.1,Tchi=.8,Te=1.,G=.001,ge=4.,
                n_p=8,n_k=8,n_mu=6,n_cm=8,n_phi=8)
    b=weak_energy_moments(M=.1,Tchi=.8,Te=1.,G=.001,ge=6.,
                n_p=8,n_k=8,n_mu=6,n_cm=8,n_phi=8)
    assert math.isclose(b["Q_chi_GeV5"]/a["Q_chi_GeV5"],
                        1.5,rel_tol=1e-14)


def test_nonfinite_couplings_fail():
    with pytest.raises(ValueError):
        same_temperature_massless_lepton_weights(float("inf"),1)
