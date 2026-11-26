"""Conservative event projection and exact FD logit interpolation tests."""
import math
import sys
from pathlib import Path
import numpy as np
import pytest
from scipy.special import expit

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from conservative_energy_projection import (
    logit_energy_interpolator,project_weak_events_energy_bins)
from relativistic_nonlinear_fermion import (
    collision_from_nonthermal_chi,distribution_from_grid)


def test_FD_exact_on_coarse_grid_and_analytic_high_energy_tail():
    M=.1;T=1.;mu=.2
    q=np.array([0.,.1,.5,2.,5.,10.,20.])
    f=expit((mu-np.hypot(q,M))/T)
    chi=logit_energy_interpolator(q,f,M)
    p=np.array([0.,.004,.25,.9,1.5,4.,8.,15.,19.,25.,33.])
    np.testing.assert_allclose(
        chi(p),expit((mu-np.hypot(p,M))/T),rtol=1e-13,atol=1e-15)


def test_PCHIP_false_equilibrium_heat_removed_on_same_very_coarse_grid():
    M=.1;T=1.
    q=np.linspace(0,32,41)
    f=expit(-np.hypot(q,M)/T)
    pchip=distribution_from_grid(q,f)
    logit_f=logit_energy_interpolator(q,f,M)
    kw=dict(M=M,Te=T,G=.001,n_p=12,n_k=12,
            n_mu=10,n_cm=12,n_phi=12,cutoff=24,p_max=24)
    bad=collision_from_nonthermal_chi(pchip,**kw)
    good=collision_from_nonthermal_chi(logit_f,**kw)
    assert abs(bad["Q_chi_GeV5"])>1e-14
    assert abs(good["Q_chi_GeV5"])<1e-17
    assert good["total_entropy_production_GeV4"]<1e-24


def test_exact_number_and_energy_projected_to_all_nodes():
    M=.1
    q=np.array([0.,.3,.7,1.4,2.2,4.,7.,11.])
    Ein=np.array([.2,.8,1.5,3.2,4.8])
    Eout=np.array([.9,2.2,.4,6.5,3.0])
    rates=np.array([-.03,.11,.2,-.08,.15])
    r=project_weak_events_energy_bins(q,M,Ein,Eout,rates)
    assert abs(r["number_residual"])<1e-14
    assert abs(r["energy_projection_defect"])<1e-13
    assert math.isclose(r["dark_energy_rate_from_grid"],
                        -r["sm_energy_rate"],abs_tol=1e-15)
    assert r["not_time_integrator"]
    assert r["not_Pauli_stable_time_update"]


def test_overflow_requires_grid_extension_and_fermi_endpoints_rejected():
    with pytest.raises(ValueError,match="outside conservative grid"):
        project_weak_events_energy_bins([0,.2,.5],.1,[.2],[3.],[.01])
    with pytest.raises(ValueError):
        logit_energy_interpolator([0,.2,.4],[.3,0.,.2],.1)
    with pytest.raises(ValueError):
        logit_energy_interpolator([0,.2,.4],[.3,1.,.2],.1)
