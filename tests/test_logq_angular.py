"""Verify endpoint-adapted angular decay quadrature independently of kinetics."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from one_to_two_fixed_bath import FlatThermalDecay
from logq_angular_quadrature import LogQDaughterKernel
from lattice_detailed_balance import EntropyBalancedLattice


def test_log_daughter_energy_change_of_variable_is_normalized():
    original=FlatThermalDecay(100,140,48)
    adapted=LogQDaughterKernel(100,140,48)
    assert adapted.angular_normalization_max_error<1e-9
    assert adapted.angular_phi_mean_energy_rel_error<1e-9
    assert np.max(abs(adapted.qphi+adapted.qchi-adapted.parentE))<1e-12
    assert adapted.min_event_daughter_q < min(original.qphi.min(),original.qchi.min())/10
    # The continuum vacuum rate does not depend on angular emission rule.
    r1=np.dot(original.rateprefactor,original.parent_fermi)
    r2=np.dot(adapted.rateprefactor,adapted.parent_fermi)
    assert abs(r1/r2-1)<1e-9


def test_logq_lattice_has_same_structural_equilibrium_and_h_theorem():
    for angular in ("mu","logq"):
        m=EntropyBalancedLattice(50,64,48,angular_rule=angular)
        out=m.checks()
        assert out["thermal_equilibrium_max_bin_rate"]<1e-11
        assert abs(out["daughter_number_rate_mismatch"])<1e-10
        assert out["free_energy_derivative"]<0
        assert out["entropy_production"]>0
        assert abs(out["dF_plus_entropy_production"])<1e-9


def test_angular_convergence_improves_at_fixed_momentum_grid():
    # This is a controlled surrogate comparison, NOT a physical abundance.
    adapted24=EntropyBalancedLattice(100,140,24,angular_rule="logq").evolve()
    adapted48=EntropyBalancedLattice(100,140,48,angular_rule="logq").evolve()
    standard24=EntropyBalancedLattice(100,140,24,angular_rule="mu").evolve()
    standard48=EntropyBalancedLattice(100,140,48,angular_rule="mu").evolve()
    delta_log=abs(adapted24["mean_phi_q"]/adapted48["mean_phi_q"]-1)
    delta_mu=abs(standard24["mean_phi_q"]/standard48["mean_phi_q"]-1)
    assert delta_log<.03
    assert delta_mu>.25
    assert delta_log<delta_mu/10
