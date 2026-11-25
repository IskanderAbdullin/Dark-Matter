"""Regression for the two complementary quartic stability criteria."""
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'kinetics'))
from direct_activation_consistency import MPL_REDUCED
from direct_pole_matching import audit
from resonant_mediator_inverse import M


def test_stationary_quartic_witness_and_finite_field_margin():
    candidate = json.loads((ROOT/'kinetics/results/direct-subthreshold-friedmann.json').read_text())[0]
    result = audit(candidate)
    stability = result['finite_field_screen']
    coupling_at_bound = result['kappa']*MPL_REDUCED
    s = 2*M**2/coupling_at_bound
    # Tightness of the sufficient quartic-only S bound can be checked
    # directly at its worst fixed-field value.
    lambda_a = stability['sufficient_lambda_A']
    a_squared = 6*(coupling_at_bound*s-M**2)/lambda_a
    positive_s = stability['lambda_S']*s**4/24
    total = positive_s + lambda_a*a_squared**2/24 + (M**2-coupling_at_bound*s)*a_squared/2
    assert abs(total/positive_s) < 2e-14
    assert np.isclose(stability['lambda_A_for_field_cap_point01_Mpl']/lambda_a, 1e-8)
    witness = result['negative_quartic_witness']
    assert witness['quartic_coefficient'] < 0.
    assert 6e5 < witness['required_over_cap'] < 8e5
    assert witness['chi_where_witness_overcomes_quadratic_over_Mpl'] > 1000
