"""Visible-electroweak-broken dark vector partial widths, without inventing Rhad."""
import math
from kinetics.phi_bose_vector_portal_visible_width_gate import (
    ME_GEV,MMU_GEV,MTAU_GEV,charged_lepton_width,
    visible_width_gate,narrow_resonant_all_charged_MB,
)
from kinetics.phi_bose_vector_portal_annihilation_resonance import (
    minimal_vector_width,narrow_width_vector_onshell_MB,
)


def test_lepton_partial_width_exact_threshold_and_high_mass_vector_limit():
    for mass in (ME_GEV,MMU_GEV,MTAU_GEV):
        assert charged_lepton_width(2*mass,.00001,mass)==0
        assert charged_lepton_width(2.01*mass,.00001,mass)>0
        rate=charged_lepton_width(100.,1e-5,mass)
        e=charged_lepton_width(100.,1e-5,ME_GEV)
        assert rate/e>.999


def test_hadron_data_are_never_fabricated_from_free_quark_thresholds():
    partial=visible_width_gate()
    assert partial["gamma_hadrons_data_driven_GeV"] is None
    assert partial["all_visible_to_electron_width_ratio"] is None
    assert partial["missing_observed_R_had_and_running_coupling"]
    lepton=partial["lepton_to_electron_width_ratio"]
    assert 2.99<lepton<3.
    specified=visible_width_gate(R_had=3.5)
    assert 6.49<specified["all_visible_to_electron_width_ratio"]<6.51
    assert specified["gamma_hadrons_data_driven_GeV"]>0
    assert (specified["dark_plus_visible_partial_GammaV_GeV"] >
            partial["dark_plus_visible_partial_GammaV_GeV"])
    assert abs(specified["total_width_relative_shift"])<1e-7


def test_resonant_NWA_branches_keep_both_dark_widths_and_redistribute_to_mu_tau():
    r=narrow_resonant_all_charged_MB()
    e=narrow_width_vector_onshell_MB(1.,.02)
    assert 2.99<r["lepton_only_relative_to_prior_eonly"]<3.001
    assert r["forward_ee_NWA_GeV4"]<e["real_intermediate_onshell_forward_GeV4"]
    assert r["forward_all_lepton_pairs_NWA_GeV4"]>2.99*r["forward_ee_NWA_GeV4"]
    assert r["forward_lepton_plus_hadron_NWA_GeV4"] is None
    assert 1e-15<r["forward_all_lepton_pairs_NWA_GeV4"]<1.1e-15
    assert minimal_vector_width()["channels_GeV"]["chi1+antichi1"]>0
    assert minimal_vector_width()["channels_GeV"]["chi2+antichi2"]>0


def test_hadronic_R_is_a_conditioned_external_input_not_a_measured_claim():
    lepton=narrow_resonant_all_charged_MB()
    hypothetical=narrow_resonant_all_charged_MB(R_had=3.5)
    assert hypothetical["forward_lepton_plus_hadron_NWA_GeV4"]>2*(
        lepton["forward_all_lepton_pairs_NWA_GeV4"])
    assert 6.49<hypothetical["visible_widths"][
        "all_visible_to_electron_width_ratio"]<6.51
    assert not hypothetical["physically_measured_R_had_used"]
    try:
        visible_width_gate(R_had=-1.)
    except ValueError:
        pass
    else: raise AssertionError("Negative physical hadronic rate accepted")
