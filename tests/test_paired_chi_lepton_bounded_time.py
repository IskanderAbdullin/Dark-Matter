"""Scoped first Pauli-bounded χ+lepton finite-volume NONLINEAR time pilot."""
import math
import numpy as np
from scipy.special import expit
from kinetics.paired_chi_lepton_bounded_time import (
    phase_space_cell_weights,quantum_entropy,system,ledger,rhs,evolve
)


def test_physical_massive_and_massless_energy_cell_integrals():
    for mass in (0.,.1,.5):
        E=np.linspace(mass,24.1+mass,21)
        w=phase_space_cell_weights(E,mass,4)
        exact=4/(6*math.pi**2)*((E[-1]*E[-1]-mass*mass)**1.5)
        assert np.all(w>0)
        # Domain starts at E=mass; compare total against exact continuum cell integral.
        assert math.isclose(sum(w),exact,rel_tol=2e-14)


def test_nonthermal_instantaneous_strong_pauli_entropy_number_and_heat():
    for mode in ("oscillatory","smooth"):
        cfg=system(bins=10,mode=mode)
        fc,fl=cfg["initial_chi"],cfg["initial_lepton"]
        dc,dl,s,err=rhs(cfg,fc,fl)
        assert s>0
        assert err<2e-25
        assert abs(np.dot(cfg["wc"],dc))<2e-24
        assert abs(np.dot(cfg["wl"],dl))<2e-24
        assert abs(np.dot(cfg["wc"],cfg["cg"]*dc)+
                   np.dot(cfg["wl"],cfg["lg"]*dl))<2e-24


def test_adaptive_fermi_bounded_time_has_rejections_entropy_and_invariants():
    s=evolve(bins=10,mode="oscillatory",
             total_time=1e10,initial_dt=1e10,atol=1e-7)
    assert s["accepted_steps"]>4
    assert s["rejected_proposals"]>=1
    assert s["min_accepted_step_entropy_gain"]>0
    assert s["entropy_gain"]>0
    assert max(abs(s["number_chi_error"]),abs(s["number_lepton_error"]),
               abs(s["total_energy_error"]))<2e-12
    assert 0<s["min_fchi"]<=s["max_fchi"]<1
    assert 0<s["min_flepton"]<=s["max_flepton"]<1
    assert s["max_step_doubling_spectral_abs_error"]<=1e-7


def test_tolerance_refinement_reduces_spectral_numerical_drift():
    a=evolve(bins=10,mode="oscillatory",atol=1e-6)
    b=evolve(bins=10,mode="oscillatory",atol=1e-7)
    c=evolve(bins=10,mode="oscillatory",atol=1e-8)
    delta_a=abs(a["fchi_first_bin"]-c["fchi_first_bin"])
    delta_b=abs(b["fchi_first_bin"]-c["fchi_first_bin"])
    assert delta_a>delta_b>0
    assert a["accepted_steps"]<b["accepted_steps"]<c["accepted_steps"]


def test_common_temperature_chemical_equilibrium_is_stationary():
    cfg=system(bins=20,mode="equal_temperature")
    dc,dl,s,err=rhs(cfg,cfg["initial_chi"],cfg["initial_lepton"])
    assert max(abs(dc))<1e-18
    assert max(abs(dl))<1e-18
    assert abs(s)<1e-30
