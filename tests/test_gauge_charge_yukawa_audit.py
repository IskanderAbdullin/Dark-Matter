"""Gauge charge selection rules for the EXACT original real-phi Yukawa."""
from kinetics.gauge_charge_yukawa_audit import analyze,finite_charge_grid


def test_real_phi_yukawa_forces_equal_vectorlike_charges():
    for i in range(-3,4):
        for j in range(-3,4):
            a=analyze(i,j,0,real_phi=True)
            assert a["renormalizable_yukawa_gauge_invariant"]==(i==j)
    assert len(finite_charge_grid())==7


def test_nonzero_chi1_gauge_charge_forces_chi2_vector_portal():
    a=analyze("1/3","1/3",0)
    assert a["z2_yukawa_even"]
    assert a["chi1_diagonal_vector_current_nonzero"]
    assert a["chi2_diagonal_vector_current_nonzero"]


def test_chi1_only_vector_current_incompatible_with_original_real_phi_yukawa():
    a=analyze(1,0,0)
    assert not a["renormalizable_yukawa_gauge_invariant"]
    assert not a["chi1_only_vector_portal_with_neutral_chi2_and_real_phi"]


def test_complex_charged_phi_restores_invariance_but_changes_the_model():
    assert analyze(1,0,-1,real_phi=False)["renormalizable_yukawa_gauge_invariant"]
    a=analyze(1,0,-1,real_phi=True)
    assert not a["one_real_phi_U1_representable"]
    assert not a["renormalizable_yukawa_gauge_invariant"]
