"""Exact count-invariant audit of restricted and extended real-scalar Yukawa."""
import numpy as np
from kinetics.phi_bose_yukawa_charge_audit import (
    SPECIES,CHANNELS,stoichiometry,analyze,
)


def test_isolated_decay_has_more_accidental_linear_counts():
    S=stoichiometry(False)
    assert S.shape==(2,5)
    assert np.linalg.matrix_rank(S)==2
    old_total=np.array([1,1,1,1,0])
    old_n1phi=np.array([1,1,0,0,1])
    u1=np.array([1,-1,1,-1,0])
    for q in (old_total,old_n1phi,u1):
        assert np.all(S@q==0)


def test_pair_channels_break_both_old_positive_dark_counts():
    S=stoichiometry(True)
    old_total=np.array([1,1,1,1,0])
    old_n1phi=np.array([1,1,0,0,1])
    assert np.array_equal(S@old_total,[0,0,-2,-2])
    assert np.array_equal(S@old_n1phi,[0,0,0,2])
    assert np.linalg.matrix_rank(S)==4
    # Exact integer 4x4 minor proves independent rank at least 4.
    assert round(np.linalg.det(S[:,[0,1,2,4]]))==-2


def test_net_dirac_U1_and_candidate_stability_Z2_survive_both_pair_channels():
    S=stoichiometry(True)
    fundamental=np.array([1,-1,1,-1,0])
    odd_Z2=np.array([1,1,0,0,1])
    assert np.all(S@fundamental==0)
    assert np.all((S@odd_Z2)%2==0)
    # Full 5-state real stoichiometric nullspace dimension exactly one,
    # so U1 spans the continuous charges for THESE four processes.
    assert len(SPECIES)-np.linalg.matrix_rank(S)==1
    assert len(CHANNELS)==4


def test_one_to_two_collision_invariants_are_not_a_full_microphysics_claim():
    a=analyze()
    assert a["decays_only_nullity"]==3
    assert a["full_nullity"]==1
    assert a["full_U1_reaction_changes"]==[0,0,0,0]
    assert a["z2_parity_changes_mod2"]==[0,0,0,0]
    assert a["old_total_fermion_reaction_changes"]==[0,0,-2,-2]
    assert a["old_N1_plus_phi_reaction_changes"]==[0,0,0,2]
