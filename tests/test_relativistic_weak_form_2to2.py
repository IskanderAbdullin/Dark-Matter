"""Independent quantitative tests of Pauli-blocked relativistic weak form."""
from __future__ import annotations
import math
import sys
from pathlib import Path
import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from relativistic_weak_form_2to2 import (
    cm_event, weak_energy_moments, rest_blocked_heat_independent)
from microphysics_recoil_heat import recoil_heating_rate_chi_at_rest


def run(n=12, cutoff=24., **kwargs):
    return weak_energy_moments(n_p=n,n_k=n,n_mu=n,n_cm=n,n_phi=n,
                               cutoff=cutoff,**kwargs)


@pytest.mark.parametrize("M,p,k,mu", [
    (.1,.3,.4,-.8),(.1,0.,.9,.2),(2.,.05,.1,.99),(.01,5.,4.,-.99),
])
def test_four_vectors_cm_azimuth_satisfy_mass_shell_and_energy(M,p,k,mu):
    u=np.array([-.91,-.5,.1,.7,.98])[:,None]
    phi=np.array([0.,.1,1.3,2.5,5.4])[None,:]
    E,Ke,t,err=cm_event(M,p,k,mu,u,phi)
    assert np.all(E>=M) and np.all(Ke>0)
    assert np.max(err)<1e-11
    assert np.all(t<=0)
    assert np.allclose(E+Ke,math.hypot(M,p)+k,rtol=2e-14,atol=1e-14)


def test_common_T_two_pauli_equilibrium_mode_by_mode_with_nonzero_mu():
    r=run(n=10,M=.1,Tchi=1.,Te=1.,muchi=.13,mue=-.2)
    assert r["Q_chi_GeV5"]==0.
    assert r["energy_strong_GeV5"]==0.
    assert r["number_strong_GeV4"]==0.
    assert all(c==0 for _,c in r["strong_modes"])
    assert r["equilibrium_control_direct_difference"]<1e-15


def test_colder_chi_heating_hotter_chi_cooling_entropy_positive():
    cold=run(n=12,M=.1,Tchi=.8,Te=1.)
    hot=run(n=12,M=.1,Tchi=1.2,Te=1.)
    assert cold["Q_chi_GeV5"]>0 and hot["Q_chi_GeV5"]<0
    for r in [cold,hot]:
        assert r["entropy_GeV4"]>0
        assert r["Q_e_GeV5"]==-r["Q_chi_GeV5"]
        assert r["max_mass_shell_defect"]<1e-11
        assert r["no_gamma_el_or_relic"]


def test_g_squared_and_internal_degeneracy_scaling():
    a=run(n=8,M=.1,Tchi=.8,Te=1.,G=.01,gchi=2,ge=2)
    b=run(n=8,M=.1,Tchi=.8,Te=1.,G=.02,gchi=2,ge=2)
    c=run(n=8,M=.1,Tchi=.8,Te=1.,G=.01,gchi=4,ge=4)
    assert math.isclose(b["Q_chi_GeV5"]/a["Q_chi_GeV5"],4,rel_tol=2e-13)
    assert math.isclose(c["Q_chi_GeV5"]/a["Q_chi_GeV5"],4,rel_tol=2e-13)


def test_independent_rest_electron_pauli_recoil_reference():
    a=rest_blocked_heat_independent(M=1.,Te=.1,G=.2,n_k=24,n_cm=24)
    b=recoil_heating_rate_chi_at_rest(1.,.1,.2)
    assert abs(a/b-1)<1e-6


def test_independent_unsymmetrised_strong_source_converges_to_weak():
    coarse=run(n=12,cutoff=24.)
    fine=run(n=24,cutoff=32.)
    assert abs(fine["weak_vs_strong_heat_relative_discrepancy"])<5e-6
    assert abs(fine["number_strong_GeV4"])<1e-9
    assert abs(fine["Q_chi_GeV5"]-coarse["Q_chi_GeV5"])/fine["Q_chi_GeV5"]<1e-3
    assert abs(coarse["number_strong_GeV4"])>abs(fine["number_strong_GeV4"])


def test_finite_cutoff_causes_real_strong_number_and_energy_leakage():
    a=run(n=24,cutoff=15.)
    b=run(n=24,cutoff=32.)
    assert abs(a["weak_vs_strong_heat_relative_discrepancy"])>1e-3
    assert abs(b["weak_vs_strong_heat_relative_discrepancy"])<5e-6
    assert abs(a["number_strong_GeV4"])>abs(b["number_strong_GeV4"])


@pytest.mark.parametrize("params", [
    {"M":0.},{"G":-1.},{"Te":0.},{"n_mu":2},
    {"Tchi":float("nan")},{"ge":0.},{"cutoff":1.},
])
def test_nonphysical_input_rejected(params):
    with pytest.raises(ValueError):
        weak_energy_moments(**params)
