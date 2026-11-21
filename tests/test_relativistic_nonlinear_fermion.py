"""Independent exact thermal controls and genuinely nonthermal Pauli source."""
from __future__ import annotations
import math
import sys
from pathlib import Path
import numpy as np
import pytest
from scipy.special import expit

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from relativistic_nonlinear_fermion import (
    collision_from_nonthermal_chi,distribution_from_grid)
from relativistic_weak_form_2to2 import weak_energy_moments


def test_arbitrary_callable_reduces_exactly_to_previous_thermal_FD_source():
    M=.1;Tchi=.8;Te=1.
    f=lambda p:expit(-np.hypot(p,M)/Tchi)
    for n in (10,16):
        a=collision_from_nonthermal_chi(
            f,M=M,Te=Te,G=.001,n_p=n,n_k=n,n_mu=n,n_cm=n,n_phi=n,
            cutoff=30.,p_max=30*Tchi)
        b=weak_energy_moments(
            M=M,Tchi=Tchi,Te=Te,G=.001,n_p=n,n_k=n,n_mu=n,n_cm=n,
            n_phi=n,cutoff=30.)
        assert math.isclose(a["Q_chi_GeV5"],b["Q_chi_GeV5"],rel_tol=1e-12)
        assert math.isclose(a["energy_strong_GeV5"],b["energy_strong_GeV5"],rel_tol=1e-10)


def test_exact_analytic_FD_at_common_T_remains_stationary_without_projection():
    M=.1;T=1.
    f=lambda p:expit((.2-np.hypot(p,M))/T)
    r=collision_from_nonthermal_chi(
        f,M=M,Te=T,mue=-.15,G=.001,n_p=12,n_k=12,n_mu=10,n_cm=12,n_phi=12)
    assert abs(r["Q_chi_GeV5"])<1e-18
    assert abs(r["number_strong_GeV4"])<1e-18
    assert max(abs(rate) for _,rate in r["strong_modes"])<1e-18
    assert r["negative_event_entropy_floor"]>=-1e-20


def make_nonth(nq=321,eps=.06,qmax=34):
    q=np.linspace(0.,qmax,nq)
    f=expit(-np.hypot(q,.1))*(1+eps*np.sin(.7*q))
    return distribution_from_grid(q,f)


def test_nonthermal_distribution_exposes_nonlinear_source_and_positive_entropy():
    r=collision_from_nonthermal_chi(
        make_nonth(),M=.1,Te=1.,G=.001,
        n_p=24,n_k=24,n_mu=16,n_cm=20,n_phi=16,cutoff=30,p_max=30)
    assert r["Q_chi_GeV5"]>4e-9
    assert r["Q_e_GeV5"]==-r["Q_chi_GeV5"]
    assert r["total_entropy_production_GeV4"]>1e-10
    assert r["negative_event_entropy_floor"]>=0
    assert abs(r["number_strong_GeV4"])<1e-10
    assert abs(r["strong_vs_weak_heat_relative"])<.001
    assert any(rate>0 for _,rate in r["strong_modes"])
    assert any(rate<0 for _,rate in r["strong_modes"])
    assert r["no_gamma_el_or_relic"]


def test_near_equilibrium_response_is_linear_Q_and_quadratic_entropy():
    a=collision_from_nonthermal_chi(
        make_nonth(eps=.03),M=.1,Te=1.,G=.001,
        n_p=16,n_k=16,n_mu=12,n_cm=16,n_phi=12,cutoff=28,p_max=28)
    b=collision_from_nonthermal_chi(
        make_nonth(eps=.06),M=.1,Te=1.,G=.001,
        n_p=16,n_k=16,n_mu=12,n_cm=16,n_phi=12,cutoff=28,p_max=28)
    assert abs(b["Q_chi_GeV5"]/a["Q_chi_GeV5"]-2)<.02
    assert abs(b["total_entropy_production_GeV4"]/a["total_entropy_production_GeV4"]-4)<.1


def test_coarse_PCHIP_fd_creates_spurious_equilibrium_heat_and_entropy():
    M=.1;T=1.
    f=lambda nq:distribution_from_grid(
        np.linspace(0.,32,nq),
        expit(-np.hypot(np.linspace(0.,32,nq),M)))
    def run(nq):
        return collision_from_nonthermal_chi(
          f(nq),M=M,Te=T,G=.001,
          n_p=20,n_k=20,n_mu=16,n_cm=18,n_phi=16,cutoff=30,p_max=30)
    coarse=run(41); fine=run(641)
    assert abs(coarse["Q_chi_GeV5"])>10*abs(fine["Q_chi_GeV5"])
    assert coarse["total_entropy_production_GeV4"]>10*fine["total_entropy_production_GeV4"]
    assert coarse["total_entropy_production_GeV4"]>0


@pytest.mark.parametrize("q,f",[
    ([0,1,1.],[.1,.2,.3]),([1,2,3],[.1,.2,.3]),
    ([0,1,2],[.1,-.2,.3]),([0,1,2],[.1,1.2,.3]),
])
def test_invalid_fermionic_tabulation_rejected(q,f):
    with pytest.raises(ValueError):
        distribution_from_grid(q,f)


def test_zero_outside_grid_and_reject_bad_callable():
    f=distribution_from_grid([0,1,2],[.5,.3,.01])
    assert float(f(3.))==0
    assert float(f(1.5))>=0
    with pytest.raises(ValueError):
        collision_from_nonthermal_chi(lambda p:-.1,n_p=3,n_k=3,n_mu=3,n_cm=3,n_phi=3)
