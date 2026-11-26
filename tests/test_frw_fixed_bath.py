"""FRW toy cosmology: comoving measure, expansion normalization, conservation."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from frw_fixed_bath_pilot import FRWFixedBathPilot
from onshell_collision_2d import OnShellCollision2D


def test_FRW_collision_prefactor_matches_fixed_bath_at_fixed_x():
    x=.1
    kappa=.02
    frw=FRWFixedBathPilot(nq=50,kappa=kappa)
    lab=OnShellCollision2D(nq=50,M1_T=x)
    y=np.zeros(100)
    dfrw=frw.rhs(x,y)
    dlab=lab.rhs(0,y)
    # d/dx = [Gamma/(x H)] d/d(t Gamma); Gamma/H=kappa*x^2.
    assert np.allclose(dfrw,dlab*(kappa*x),rtol=1e-11,atol=1e-12)


def test_FRW_exact_daughter_invariant_at_varying_temperature():
    frw=FRWFixedBathPilot(nq=50,kappa=.03)
    for x in (.1,.5,1.,3.,10.):
        y=np.r_[np.exp(-frw.e),0.4*np.exp(-frw.e)]
        d=frw.rhs(x,y)
        mismatch=np.dot(frw.vphi,d[:frw.nq])-np.dot(frw.vchi,d[frw.nq:])
        assert abs(mismatch)<1e-10


@pytest.mark.parametrize("x",[.1,1.,5.])
def test_FRW_pairwise_Bose_Fermi_detailed_balance(x):
    frw=FRWFixedBathPilot(nq=50,kappa=.01)
    y=np.r_[1/np.expm1(frw.e+.3),1/(1+np.exp(frw.e-.3))]
    d=frw.rhs(x,y)
    assert np.max(abs(d[:frw.nq]*frw.vphi))<1e-11
    assert np.max(abs(d[frw.nq:]*frw.vchi))<1e-11


def test_radiation_dominated_source_suppressed_at_late_x():
    frw=FRWFixedBathPilot(nq=40,kappa=.001)
    result=frw.integrate(x_start=.05,x_mid=20,x_end=30)
    assert result["integration"]["success"]
    assert result["frozen_number_to_1e-4_over_late_window"]
    assert abs(result["final"]["number_difference_over_T3"])<2e-7
    # Weak-coupling toy produces a hot, not instantaneous-burst-cold spectrum.
    assert result["final"]["mean_phi_q_over_T"]>1.


def test_coupling_scan_changes_frozen_spectrum_not_just_amplitude():
    weak=FRWFixedBathPilot(nq=40,kappa=.001).integrate()["final"]
    strong=FRWFixedBathPilot(nq=40,kappa=.1).integrate()["final"]
    assert strong["nphi_over_T3"]>weak["nphi_over_T3"]
    assert strong["mean_phi_q_over_T"]<weak["mean_phi_q_over_T"]
