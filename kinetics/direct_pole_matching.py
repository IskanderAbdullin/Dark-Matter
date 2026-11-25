"""Stability and finite matching tests for the even, renormalizable A-S-chi potential."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid

from direct_activation_consistency import MPL_REDUCED
from early_bose_bootstrap import N0, rearrange
from late_activation_overlap import H0, Or
from direct_subthreshold_activation import SCALE, sunset_threshold
from quartic_consistency import lambda_from_sigma, sigma_over_m
from resonant_mediator_inverse import HBAR, M
from soft_mixing_activation import variance_table


def audit(candidate):
    p = candidate['params']
    v = p['vev_eV']
    mu2 = p['mu_eV']**2
    chi_quartic = p['lambda_chi']
    detuning = candidate['a']*SCALE
    ms2 = 4*M*M-detuning
    g2 = candidate['g2_eV2']
    kappa = np.sqrt(g2)/v

    # Mixed convention: V contains h_ij phi_i^2 phi_j^2/4, and
    # h_AS radiatively generates A^4 at Q=mS/2 relative to Q=mS.
    # Even this very generous A cap has failed collisionless evolution.
    a_quartic_cap = float(lambda_from_sigma(1e-7))
    h_as_cap = np.sqrt(a_quartic_cap*16*np.pi**2/(3*np.log(2.)))
    h_a_chi_cap = 32*np.pi**2*mu2/M**2
    h_s_chi_cap = 32*np.pi**2*mu2/ms2

    # Certified negative quartic witness, taking every positive coefficient
    # at its generous cap. Choose S/chi near the maximum of h^2/(a+b*h^4).
    s_quartic_max = 4*np.pi
    h = (6*chi_quartic/s_quartic_max)**.25
    negative_a2 = kappa*h/2 - (h_a_chi_cap + h_as_cap*h*h)/4
    positive_constant = chi_quartic/4 + h_s_chi_cap*h*h/4+s_quartic_max*h**4/24
    min_a_squared = 12*negative_a2/a_quartic_cap
    quartic_at_witness = positive_constant - 6*negative_a2**2/a_quartic_cap
    quadratic_at_witness = (M*M*min_a_squared + ms2*h*h-mu2)/2
    a_quartic_needed = 6*negative_a2**2/positive_constant

    # Positive masses allow a much weaker, constructive finite-|chi| test.
    # Minimize V over A^2 and use V_S >= lambda_S S^4/24.
    # For g=|kappa chi|, max_{S>M^2/g} (g S-M^2)^2/S^4 = g^4/(16 M^4).
    s_threshold_budget = .1
    s_quartic = s_threshold_budget*32*np.pi**2*detuning/ms2
    finite_field_a = 9*(kappa*MPL_REDUCED)**4/(4*s_quartic*M**4)
    constant_a = 3*g2*candidate['full_friedmann']['max_field']**2/ms2
    finite_thresholds = []
    for factor in (.5, 1., 2.):
        terms = sunset_threshold(g2, factor*M)
        # With a constant A^4 term the A figure-8 is chi independent.
        new = terms['total_eV4']-terms['figure8_eV4']
        finite_thresholds.append(dict(Q_over_M=factor,
                                      original_ratio=2*abs(terms['total_eV4'])/p['X'],
                                      constant_A4_ratio=2*abs(new)/p['X']))
    return dict(mesh=candidate['mesh'], kappa=float(kappa), detuning_eV2=detuning,
                naturalness_test_caps=dict(lambda_A=a_quartic_cap, h_AS=float(h_as_cap),
                                           h_A_chi=float(h_a_chi_cap), h_S_chi=float(h_s_chi_cap),
                                           lambda_S=float(s_quartic_max)),
                negative_quartic_witness=dict(S_over_chi=float(h),
                    A_squared_over_chi_squared=float(min_a_squared),
                    quartic_coefficient=float(quartic_at_witness),
                    necessary_lambda_A_at_fixed_portal_caps=float(a_quartic_needed),
                    required_over_cap=float(a_quartic_needed/a_quartic_cap),
                    chi_where_witness_overcomes_quadratic_over_Mpl=float(
                        np.sqrt(quadratic_at_witness/(-quartic_at_witness))/MPL_REDUCED)),
                finite_field_screen=dict(chi_cap_over_Mpl=1., lambda_S=float(s_quartic),
                    one_loop_S_mass_threshold_over_detuning=s_threshold_budget,
                    lambda_A_for_field_cap_point01_Mpl=float(finite_field_a*1e-8),
                    sufficient_lambda_A=float(finite_field_a),
                    sufficient_sigma_cm2_g=float(sigma_over_m(finite_field_a)),
                    original_trajectory_square_equivalent_constant_A4=float(constant_a),
                    trajectory_square_equivalent_sigma_cm2_g=float(sigma_over_m(constant_a)),
                    finite_field_range_over_Mpl_with_trajectory_A4=float((constant_a/finite_field_a)**.25),
                    finite_trigger_thresholds=finite_thresholds),
                scope='Positive even portals, selected threshold budgets without cancellations, and two-field boundedness at fixed chi. No UV completion, global quantum stability, cosmological production or abundance proof.')


def elastic_bound(row, coupling):
    """Conservative contact-only gain bound with prescribed A caps.

    Assume H >= 0.5 H0 sqrt(Or) r^2, leaving margin for SM entropy release.
    The finite-amplitude tree contact result does not cover a negative mass
    shift, coherence, independent species or a cap violation by the portal.
    """
    f = np.asarray(row['frozen_cell_spectrum']['f_phi'])
    occupation_cap = 2.1*float(np.max(f))
    number_prefactor = 2*N0
    a = number_prefactor/M
    pstar_today = np.cbrt(6*np.pi**2*number_prefactor/occupation_cap)
    b = occupation_cap*pstar_today**2/(4*np.pi**2)
    h = .5*H0*np.sqrt(Or)*HBAR
    growth = 1+2*occupation_cap
    # Integral over all 0 < r < infinity of min(a^2*r, b^2/r)
    # with measure d ln r equals 2*a*b. Divide by the reference comoving N.
    l1_bound = 4*coupling**2*growth*2*a*b/(32*np.pi*N0*h)
    _, temperature = variance_table(row)
    start = row['parent_GeV']*1e9/(.1*temperature)
    r = np.geomspace(1., start, 4096)
    intensity, _ = rearrange(number_prefactor*r**3, occupation_cap, M)
    gain = coupling**2*growth*intensity**2/(64*np.pi**3*M*h*r*r)
    p_peak_today = row['final']['peak_bin_q']*temperature
    p_min_today = row['frozen_cell_spectrum']['edges'][0]*temperature
    return dict(occupation_cap=occupation_cap, number_cap_over_reference=2.,
                radiation_H_lower_bound_coefficient=.5,
                initial_redshift_approx=float(start),
                lambda_times_occupation_cap=float(coupling*occupation_cap),
                comoving_L1_bound_over_reference=float(l1_bound),
                positive_integrated_gain_per_mode=float(trapezoid(gain, np.log(r))),
                hartree_squared_gap_over_peak_momentum_squared_upper=float(coupling*b/(2*p_peak_today**2)),
                hartree_squared_gap_over_minimum_momentum_squared_upper=float(coupling*b/(2*p_min_today**2)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source = Path(__file__).with_name('results')/'direct-subthreshold-friedmann.json'
    production = json.loads((source.parent/'converged-independent-kernels-2026-11-25.json').read_text(encoding='utf-8-sig'))
    results = []
    for candidate in json.loads(source.read_text()):
        report = audit(candidate)
        row = production[f'parent_angle{candidate["mesh"]}_sparse_k00685'][0]
        field = report['finite_field_screen']
        report['constant_contact_bounds'] = {
            'small_field_range': elastic_bound(row, field['lambda_A_for_field_cap_point01_Mpl']),
            'original_square_equivalent': elastic_bound(row, field['original_trajectory_square_equivalent_constant_A4']),
            'planck_field_range': elastic_bound(row, field['sufficient_lambda_A']),
        }
        results.append(report)
    args.output.write_text(json.dumps(results, indent=2, allow_nan=False)+'\n')
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
