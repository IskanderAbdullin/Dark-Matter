"""Finite-mass daughter: exact two-body laboratory source, energy and IR limit."""
from __future__ import annotations

import math
import sys
from pathlib import Path
import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import expit

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from massive_daughter_decay_source import (
    rest_kinematics, massive_daughter_vacuum_spectrum,
    integrated_source_and_energy, zero_momentum_bose_occupation_source, run,
)
from onshell_collision_2d import vacuum_reference
from best_vacuum_spectrum import vacuum_phi_source_differential


@pytest.mark.parametrize("q",[1e-4,.003,.05,.1,.5,2.])
def test_zero_mass_limit_recovers_independent_previous_exact_spectral_formula(q):
    new=massive_daughter_vacuum_spectrum(q,mphi=0.,mchi=0.)
    old=vacuum_phi_source_differential(q)
    assert math.isclose(new,old,rel_tol=1e-11,abs_tol=1e-15)


@pytest.mark.parametrize("mphi,mchi,q",[
    (.001,0,.001),(.001,0,.1),(.01,0,.03),
    (.02,.015,.08),(.04,.02,.5)])
def test_exact_massive_spectrum_matches_independently_integrated_parent_FD_energy(
        mphi,mchi,q):
    M=.1
    k=rest_kinematics(M,mphi,mchi)
    estar=k["rest_phi_energy"];pstar=k["rest_daughter_momentum"]
    lab=math.hypot(q,mphi)
    Emin=M*(estar**2+q**2)/(lab*estar+q*pstar)
    Emax=M*(lab*estar+q*pstar)/(mphi*mphi)
    direct=quad(lambda E:expit(-E),Emin,Emax,epsabs=1e-12,epsrel=1e-12)[0]
    direct*=4*M*M/(4*math.pi**2*pstar)*(q/lab)
    exact=massive_daughter_vacuum_spectrum(q,M,mphi,mchi)
    assert math.isclose(exact,direct,rel_tol=1e-10,abs_tol=1e-13)


@pytest.mark.parametrize("mass_ratio",[0.,.001,.005,.01,.02,.05,.1,.2])
def test_fixed_parent_width_makes_total_number_independent_of_daughter_masses(mass_ratio):
    out=integrated_source_and_energy(mphi=.1*mass_ratio)
    assert abs(out["number_integral_relative_error"])<1e-9
    assert abs(out["phi_physical_energy_relative_error"])<1e-9
    assert math.isclose(out["total_phi_number_rate"],vacuum_reference(),rel_tol=1e-9)
    assert out["no_Bose_stimulation_or_zero_mode_evolution"]


def test_finite_mass_changes_energy_partition_even_if_total_rate_stays_constant():
    a=integrated_source_and_energy(mphi=0.)
    b=integrated_source_and_energy(mphi=.02)
    assert math.isclose(a["kinematics"]["phi_energy_fraction"],.5,rel_tol=1e-14)
    assert math.isclose(b["kinematics"]["phi_energy_fraction"],.52,rel_tol=1e-14)
    assert math.isclose(b["kinematics"]["chi_energy_fraction"],.48,rel_tol=1e-14)
    assert math.isclose(a["total_phi_number_rate"],b["total_phi_number_rate"],rel_tol=1e-9)


@pytest.mark.parametrize("mphi",[.0005,.001,.002,.01,.02])
def test_finite_mass_ir_mode_source_has_exact_nonzero_limit(mphi):
    q=1e-7
    value=massive_daughter_vacuum_spectrum(q,mphi=mphi)
    occupation_rate=2*math.pi**2*value/q**2
    exact=zero_momentum_bose_occupation_source(mphi=mphi)
    assert exact>0
    assert math.isclose(occupation_rate,exact,rel_tol=5e-6)


def test_ir_limit_and_vacuum_number_density_are_nonuniform_in_massless_limit():
    assert massive_daughter_vacuum_spectrum(0.,mphi=.001)==0.
    assert massive_daughter_vacuum_spectrum(0.,mphi=0.)==0.
    assert zero_momentum_bose_occupation_source(mphi=0.)==0.
    c=zero_momentum_bose_occupation_source(mphi=.001)
    assert 26000<c<28000
    assert zero_momentum_bose_occupation_source(mphi=.1*.0001)<1e-100


@pytest.mark.parametrize("bad",[
    {"M":0.},{"mphi":-.001},{"mphi":.1},
    {"mphi":.07,"mchi":.04},{"mchi":.1},{"M":float("inf")},])
def test_closed_threshold_and_unphysical_masses_rejected(bad):
    with pytest.raises(ValueError):
        rest_kinematics(**bad)


def test_invalid_momenta_or_negative_width_rejected():
    with pytest.raises(ValueError):
        massive_daughter_vacuum_spectrum([-1,.01])
    with pytest.raises(ValueError):
        massive_daughter_vacuum_spectrum(.1,rest_width=0)
    with pytest.raises(ValueError):
        zero_momentum_bose_occupation_source(mphi=.001,gphi=0)


def test_massive_source_survey_is_not_misclassified_as_primordial_DM():
    out=run()
    assert len(out["cases"])==8
    assert "No nonlinear condensate" in out["research_scope"]
    assert out["cases"][0]["finite_mass_q0_Bose_occupation_source"]==0.
    assert max(x["finite_mass_q0_Bose_occupation_source"]
               for x in out["cases"])>10000
