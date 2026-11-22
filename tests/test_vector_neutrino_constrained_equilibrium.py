"""Constrained entropy maximum of topologically corrected Zprime-neutrino graph."""
import math
import numpy as np
import pytest
from kinetics.vector_neutrino_shell_radial_endpoints import adjusted_graph
from kinetics.vector_neutrino_reversible_event_graph import (
  synthetic_spectra,ledger,sweep,event_moments,
)
from kinetics.vector_neutrino_constrained_equilibrium import (
  constrained_equilibrium,long_time_pilot,
)


@pytest.mark.parametrize("h",(.4,.2,.1))
def test_independent_three_charge_FD_BE_entropy_maximum(h):
    g=adjusted_graph(.4,h,4,scheme="shell_cell_integral")
    initial=synthetic_spectra(g)
    info,(b,f,a)=constrained_equilibrium(g,initial)
    assert info["residual_max"]<1e-12
    assert info["entropy_gain_max"]>0
    assert abs(info["mu_nu_over_T"]-info["mu_antinu_over_T"])>.01
    assert info["nonzero_lepton_asymmetry"]
    before=ledger(g,b,f,a)
    B,F,A=sweep(g,b,f,a,dt=3)
    assert max(np.max(abs(B-b)),np.max(abs(F-f)),np.max(abs(A-a)))<3e-13
    after=ledger(g,B,F,A)
    assert abs(after["entropy"]-before["entropy"])<2e-12
    assert abs(event_moments(g,b,f,a)["dNV"])<1e-12


def test_long_time_distinct_neutrino_modes_approach_independent_max_entropy():
    r=long_time_pilot(h=.2,tau=500,dt=1)
    assert r["microreversible_edges"]>=200
    assert r["max_occupation_abs_error_to_FD_BE_equilibrium"]<.002
    assert 0<r["entropy_gap_to_independent_constrained_max"]<3e-6
    assert r["entropy_gain"]>0 and r["min_step_entropy_gain"]>0
    for key in ("conserved_C_residual","conserved_L_residual",
                "conserved_E_residual"):
        assert abs(r[key])<2e-12


def test_first_order_time_bias_persists_despite_exact_invariants():
    r1=long_time_pilot(h=.4,tau=100,dt=2)
    r2=long_time_pilot(h=.4,tau=100,dt=.25)
    assert abs(r1["lowest_vector_b"]-r2["lowest_vector_b"])>.005
    assert r2["entropy_gap_to_independent_constrained_max"]>0
    assert r2["max_occupation_abs_error_to_FD_BE_equilibrium"]>0


def test_reject_incommensurate_timestep():
    with pytest.raises(ValueError):
        long_time_pilot(h=.4,tau=1,dt=.3)
