"""Independent two-energy integration tests for quantum 1<->2 chemical affinity."""
from __future__ import annotations

import math
import sys
from pathlib import Path
import numpy as np
import pytest
from scipy.integrate import quad

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from thermal_chemical_affinity import (
    thermal_1to2_spectrum, integrated_thermal_source, independent_2d_integrand,
    independent_physical_energy_reference, run,
)
from onshell_collision_2d import vacuum_reference


@pytest.mark.parametrize("q",[.003,.02,.1,1.])
def test_analytic_finite_bose_fermi_source_matches_independent_2d_gain_loss(q):
    m=.1
    ymin=m*m/(4*q)
    numerical=quad(
        lambda y: independent_2d_integrand(q,y,mu1=-.15,mu2=.1,mu_phi=-.3),
        ymin,np.inf,epsabs=1e-13,epsrel=1e-12,limit=300)[0]
    numerical*=4*m/(2*math.pi**2)
    exact=thermal_1to2_spectrum(q)
    assert math.isclose(exact,numerical,rel_tol=1e-10,abs_tol=1e-13)


def test_nonzero_bose_population_inverse_decay_and_chemical_affinity_sign():
    positive=integrated_thermal_source(mu1=-.15)
    negative=integrated_thermal_source(mu1=-.25)
    equilibrium=integrated_thermal_source(mu1=-.2)
    assert math.isclose(positive["n_phi_net_rate"],.0008565027253502168,rel_tol=1e-9)
    assert math.isclose(negative["n_phi_net_rate"],-.000826365040388528,rel_tol=1e-9)
    assert abs(equilibrium["n_phi_net_rate"])<1e-14
    assert positive["n_phi_net_rate"]>0
    assert negative["n_phi_net_rate"]<0
    assert positive["net_rate_sign_matches_affinity"]
    assert negative["net_rate_sign_matches_affinity"]
    assert abs(positive["n_phi_net_rate"]+positive["n_chi1_net_rate"])<1e-16
    assert positive["n_phi_net_rate"]==positive["n_chi2_net_rate"]


def test_exact_thermal_equilibrium_spectrum_zero_pointwise():
    q=np.geomspace(1e-5,30.,250)
    spectrum=thermal_1to2_spectrum(q,mu1=.1-.3)
    assert np.max(abs(spectrum))<1e-15


def test_vacuum_limit_recovers_old_independent_one_dimensional_reference():
    near_vacuum=integrated_thermal_source(mu1=0.,mu2=-30.,mu_phi=-30.)
    assert math.isclose(near_vacuum["n_phi_net_rate"],vacuum_reference(),
                        rel_tol=1e-9)
    assert math.isclose(near_vacuum["n_phi_net_rate"],.0164959084596,rel_tol=1e-9)


def test_equal_A_B_analytic_limit_avoids_numerical_zero_division():
    # Exactly A=B at q=.5, mu1=.6, mu2=.1.
    q=.5
    mu1=.6
    ymin=.1**2/(4*q)
    direct=quad(lambda y:independent_2d_integrand(
        q,y,mu1=mu1,mu2=.1,mu_phi=-.3),ymin,np.inf,
        epsabs=1e-12,epsrel=1e-12)[0]
    predicted=thermal_1to2_spectrum(q,mu1=mu1)
    assert math.isclose(predicted,4*.1/(2*math.pi**2)*direct,
                        rel_tol=1e-10)


def test_vectorization_and_ir_boundary_remain_finite_nonnegative_for_positive_affinity():
    q=np.r_[0.,np.geomspace(1e-6,20,100)]
    vals=thermal_1to2_spectrum(q)
    assert vals.shape==q.shape
    assert vals[0]==0
    assert np.all(np.isfinite(vals))
    assert np.all(vals>=0)
    assert thermal_1to2_spectrum(0.)==0.


@pytest.mark.parametrize("kwargs",[
    {"mass":0.},{"g1":0},{"rest_width":0.},{"mu_phi":0},
    {"q":-1.},{"q":np.nan}])
def test_invalid_kinematics_or_bose_chemical_potential_fail_closed(kwargs):
    kwargs={"q":.1,**kwargs}
    with pytest.raises(ValueError):
        thermal_1to2_spectrum(**kwargs)


def test_cli_report_not_misrepresented_as_BEST_or_cosmological_freezeout():
    record=run()
    assert not record["external_BEST_executed"]
    assert len(record["cases"])==3
    assert record["vacuum_relative_error"]<1e-9
    assert "no full nonlinear time evolution" in record["warning"]



def test_all_three_physical_energy_fluxes_from_original_2d_quantum_integrand():
    r=independent_physical_energy_reference()
    d=r["physical_energy_derivatives"]
    assert math.isclose(d["phi"],.0006730451878862239,rel_tol=1e-7)
    assert math.isclose(d["chi2"],.0011371530981677406,rel_tol=1e-7)
    assert math.isclose(d["chi1"],-.0018101982860539644,rel_tol=1e-7)
    assert abs(r["energy_conservation_residual"])<1e-8
    assert abs(r["analytic_vs_independent_2d_phi_energy_error"])<1e-8
    assert r["two_dimensional_original_quantum_factors_used"]


def test_inverse_dominated_chemical_snapshot_reverses_all_energy_flows():
    forward=independent_physical_energy_reference(mu1=-.15)
    backward=independent_physical_energy_reference(mu1=-.25)
    for k in ("chi1","chi2","phi"):
        assert forward["physical_energy_derivatives"][k]*backward["physical_energy_derivatives"][k]<0
    assert abs(backward["energy_conservation_residual"])<1e-8


def test_energy_reference_rejects_chemical_potentials_too_big_for_cutoff():
    with pytest.raises(ValueError,match="requires"):
        independent_physical_energy_reference(mu1=10.)
