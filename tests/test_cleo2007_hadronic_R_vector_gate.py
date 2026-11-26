"""Direct, error-resolved CLEO III 2007 continuum hadronic R(s) for A'10."""
import math
from kinetics.cleo2007_hadronic_R_reference import (
    CLEO2007_DATA, measured_R_at_or_near,
)
from kinetics.phi_bose_vector_portal_visible_width_gate import (
    CLEO2007_measured_visible_gate,narrow_resonant_all_charged_MB,
)


def test_table_VII_seven_primary_experimental_rows_and_uncertainty_correlations():
    assert len(CLEO2007_DATA)==7
    x=measured_R_at_or_near(10.)
    assert x["physical_energy_GeV"]==9.996
    assert x["R_had_continuum_Born"]==3.497
    assert x["stat_error"]==.004
    assert x["correlated_systematic_error"]==.064
    assert x["uncorrelated_systematic_error"]==.043
    assert math.isclose(x["total_one_point_quadrature_error"],
                        math.sqrt(.004**2+.064**2+.043**2),rel_tol=1e-14)
    assert math.isclose(x["mass_gap_GeV"],.004,abs_tol=1e-14)
    assert x["experimental_R_at_nearby_not_exact_mass"]
    assert x["relative_one_point_error"]>.02
    assert "10.1103/PhysRevD.76.072008"==x["DOI"]


def test_offpoint_mass_is_explicit_physical_extrapolation_not_fake_exact_R():
    try:measured_R_at_or_near(10.,max_gap_GeV=.0039)
    except ValueError:pass
    else:raise AssertionError("4MeV mismatch hidden")
    assert measured_R_at_or_near(9.996,0.)[
        "experimental_R_at_nearby_not_exact_mass"] is False
    try:measured_R_at_or_near(11.)
    except ValueError:pass
    else:raise AssertionError("Far outside experiment accepted")


def test_measured_R_propagates_to_dark_vector_partial_SM_width_and_NWA():
    a=CLEO2007_measured_visible_gate()
    d=a["measured_primary"]
    assert d["R_had_continuum_Born"]==3.497
    lo,hi=a["measured_onepoint_width_range_GeV"]
    assert lo<a["partial_width_GeV"]<hi
    qlo,qhi=a["measured_onepoint_NWA_range_GeV4"]
    assert qlo<a["NWA_chi2_pair_to_all_visible_T1_forward_GeV4"]<qhi
    assert qlo>0
    assert 6.4<a["ratio_NWA_all_visible_to_prior_eonly"]<6.6
    assert a["nonexperimental_4MeV_mediator_interpolation_uncertainty_omitted"]
    print("CLEO2007_DATA_DRIVEN_A_PRIME_WIDTH",
          repr(dict(experimental_R=d["R_had_continuum_Born"],
            R_error=d["total_one_point_quadrature_error"],
            measured_energy=d["physical_energy_GeV"],
            partial_width=a["partial_width_GeV"],
            partial_width_range=a["measured_onepoint_width_range_GeV"],
            all_visible_NWA_gamma=a[
               "NWA_chi2_pair_to_all_visible_T1_forward_GeV4"],
            all_visible_NWA_range=a[
               "measured_onepoint_NWA_range_GeV4"],
            ratio_to_previous_eonly=a["ratio_NWA_all_visible_to_prior_eonly"])),
          flush=True)


def test_nearby_mediator_energy_not_physical_cross_section_uniformity_claim():
    x=measured_R_at_or_near(10.)
    assert x["mass_gap_GeV"]>0
    assert x["continuum_threshold_resonance_and_vacuum_only"]
    assert x["no_dark_vector_model_fit"]
