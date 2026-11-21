"""Vacuum QED exact σv controls; FD e± initial LOSS only, not full operator."""
import math
import pytest
from kinetics.qed_pair_annihilation_control import (
    ALPHA,sigma_v_cm,sigma_v_cm_angular,sigma_v_geant_rest_reference,
    initial_fd_loss_density,initial_fd_density,quantum_event_balance)

@pytest.mark.parametrize("beta",[.01,.1,.5,.8,.95])
def test_qed_exact_sigma_vs_independent_amplitude_and_geant(beta):
    s=4/(1-beta**2)
    closed=float(sigma_v_cm(s))
    angular=sigma_v_cm_angular(s,n=160)
    heitler=sigma_v_geant_rest_reference(s)
    assert math.isclose(closed,angular,rel_tol=1e-11)
    assert math.isclose(closed,heitler,rel_tol=1e-11)

@pytest.mark.parametrize("beta",[0.,1e-6,1e-4])
def test_vacuum_threshold_and_identical_photon_count(beta):
    s=4/(1-beta**2)
    assert math.isclose(float(sigma_v_cm(s)),math.pi*ALPHA**2,
                        rel_tol=2e-8)

def test_low_temperature_initial_loss_is_nr_consistent_but_not_net():
    T=.025
    d=initial_fd_density(T)
    n16=initial_fd_loss_density(T,n=16)
    n32=initial_fd_loss_density(T,n=32)
    ratio=n32/d**2/(math.pi*ALPHA**2)
    assert math.isclose(n16,n32,rel_tol=1e-8)
    assert .98<ratio<1.
    assert ratio<1. # nonzero thermal corrections
    # FD loss alone DOES NOT represent equilibrium net zero; reverse omitted.

def test_local_quantum_balance_chemical_affinity_and_entropy():
    eq=quantum_event_balance(2,3,2.5,2.5,1,.2,-.2)
    assert math.isclose(eq["forward"],eq["reverse"],rel_tol=2e-14)
    assert abs(eq["affinity_over_T"])<1e-15
    shifted=quantum_event_balance(2,3,2.5,2.5,1,.21,-.2)
    assert math.isclose(shifted["ratio"],math.exp(.01),rel_tol=1e-12)
    assert shifted["entropy_event_factor"]>0
