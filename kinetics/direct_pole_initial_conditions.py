"""Conditional trigger-seed and coherent mediator misalignment diagnostics."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import ellipk, gamma

from coherent_floquet import oscillator, peak
from direct_activation_consistency import CGS_RHO_TO_EV4, MPL_REDUCED, T0
from direct_subthreshold_friedmann import solve_friedmann
from quartic_consistency import MPL
from resonant_mediator_inverse import HBAR, M


def rd_trajectory(mass, quartic, amplitude, match_time, *, rtol=1e-9):
    """Exact homogeneous KG through onset for constant g* and radiation H=1/(2t)."""
    omega0 = np.sqrt(mass**2+quartic*amplitude**2/6)
    mu2 = (mass/omega0)**2
    t0 = 1e-6

    def rhs(t, state):
        x, v = state
        return v, -1.5*v/t - mu2*x - (1-mu2)*x**3

    sol = solve_ivp(rhs, (t0, match_time), [1-t0*t0/5, -2*t0/5],
                    rtol=rtol, atol=rtol/20, method='DOP853',
                    dense_output=True)
    if not sol.success:
        raise RuntimeError(sol.message)
    return omega0, sol


def mediator_relic(amplitude, mass, quartic, match_time=800., rtol=1e-9):
    omega0, sol = rd_trajectory(mass, quartic, amplitude, match_time, rtol=rtol)
    times = np.linspace(.9*match_time, match_time, 40)
    x, v = sol.sol(times)
    fields = amplitude*x
    energies = (omega0*amplitude*v)**2/2 + mass**2*fields**2/2 + quartic*fields**4/24
    actions = [oscillator(energy, mass, quartic) for energy in energies]
    hubble = omega0/(2*times)
    temperature = np.sqrt(hubble*MPL*1e9/(1.66*np.sqrt(106.75)))
    redshift = temperature/T0*(106.75/3.909)**(1/3)
    rho_h2 = 1.8783416e-29*CGS_RHO_TO_EV4
    estimates = np.array([a['action'] for a in actions])*mass/redshift**3/rho_h2
    last = actions[-1]
    transition_amplitude = np.sqrt(12)*mass/np.sqrt(quartic)
    approximate_quartic_cycles = (
        .847213*last['frequency_scale']/hubble[-1]
        * max(np.sqrt(last['Phi2'])/transition_amplitude-1, 0.)/(2*np.pi)
    )
    return dict(amplitude_eV=amplitude, matching_time=match_time,
                omega_h2_coherent=float(np.median(estimates)),
                phase_span_relative=float(np.ptp(estimates)/np.median(estimates)),
                H_match_eV=float(hubble[-1]), temperature_match_GeV=float(temperature[-1]/1e9),
                quartic_cycles_after_match_before_crossover=float(approximate_quartic_cycles),
                elliptic_parameter_at_match=float(last['elliptic_parameter']))


def conformal_quartic_reference(mass, quartic, amplitude):
    # In radiation domination a'=constant. For an initially frozen S, the
    # conformal field a*S starts at zero with derivative a'*amplitude.
    # Its pure-quartic conformal energy is therefore (a'*amplitude)^2/2.
    a_prime = 1.66*np.sqrt(106.75)*T0**2*(3.909/106.75)**(2/3)/(MPL*1e9)
    action = oscillator((a_prime*amplitude)**2/2, 0., quartic)['action']
    return mass*action/(1.8783416e-29*CGS_RHO_TO_EV4)


def spectator_scan(mass, quartic, quartic_floquet):
    finite_random_walk = np.sqrt(60)/(2*np.pi)
    equilibrium = (9/(np.pi**2*quartic))**.25*np.sqrt(gamma(.75)/gamma(.25))
    regimes = []
    for budget in (.0012, .012, .12):
        ref = 1e22*(budget/conformal_quartic_reference(mass, quartic, 1e22))**(2/3)
        amp = np.exp(brentq(lambda l: mediator_relic(np.exp(l), mass, quartic)['omega_h2_coherent']/budget-1,
                            np.log(ref/3), np.log(ref*3), xtol=1e-7))
        measured = mediator_relic(amp, mass, quartic)
        checked = mediator_relic(amp, mass, quartic, 2500., rtol=2e-10)
        regimes.append(dict(budget_omega_h2=budget, **measured,
                            fine_match=checked,
                            omega_h2_independent_conformal_reference=float(conformal_quartic_reference(mass, quartic, amp)),
                            pure_quartic_peak_amplitude_efolds_after_matching_estimate=float(
                                measured['quartic_cycles_after_match_before_crossover']
                                *4*ellipk(.5)*quartic_floquet['mu']),
                            H_inf_GeV_for_60_efold_rms=float(amp/finite_random_walk/1e9),
                            H_inf_GeV_for_equilibrium_rms=float(amp/equilibrium/1e9)))
    return dict(mass_eV=mass, quartic=quartic, gstar_matching=106.75,
                H_inf_rms_factors=dict(random_walk_60=float(finite_random_walk),
                                      quartic_stationary=float(equilibrium)),
                pure_quartic_linear_floquet=quartic_floquet, budgets=regimes,
                scope='Instant radiation domination, constant high-T SM degrees at matching, homogeneous coherent branch only. Self-resonance and fragmentation can invalidate these otherwise adiabatic numbers.')


def trigger_scan(candidate, row):
    p = candidate['params']
    seed = p['seed']
    out = []
    for fac in (.1, .5, np.exp(-.01), 1., np.exp(.01), 2., 10.):
        q = dict(p, seed=seed*fac)
        summary, history = solve_friedmann(row, q, rtol=2e-9, n=7000,
                                          order=32, return_history=True)
        hit = np.flatnonzero(np.abs(history['x'])**4 >= .5)[0]
        span = slice(hit-1, hit+1)
        ln_a = np.interp(.5, np.abs(history['x'][span])**4, history['log_a'][span])
        h = np.interp(ln_a, history['log_a'][span], history['H'][span])*HBAR
        out.append(dict(seed_over_reference=float(fac), first_z_half=float(np.expm1(-ln_a)),
                        max_extra_fraction=summary['max_extra_fraction_zle100'],
                        H0_ratio=summary['H0_ratio'], H_at_half_eV=float(h),
                        tau=float(candidate['sigma0']*summary['v0']*summary['cosmic']
                                  * 0.12*1.8783416e-29)))
    response = (out[4]['first_z_half']-out[2]['first_z_half'])/.02
    v = p['vev_eV']
    wall_tension = 2*np.sqrt(2)*p['mu_eV']*v*v/3
    return dict(seed_chi_eV=seed*v, seed_trials=out, dz_half_d_log_seed=response,
                H_inf_GeV_for_10_percent_seed_fluctuation_single_mode=float(.1*seed*v*2*np.pi/1e9),
                H_inf_GeV_for_10_percent_seed_fluctuation_60_efold_rms=float(.1*seed*v*2*np.pi/np.sqrt(60)/1e9),
                wall_tension_eV3_unmodified_double_well=float(wall_tension),
                one_wall_per_hubble_fraction_at_half=float(
                    wall_tension/(3*MPL_REDUCED**2*out[3]['H_at_half_eV'])),
                eta_chi_quartic_strength=float((p['m_eta_at_vev_eV']/v)**2),
                scope='Separate homogeneous patches with fixed overdamped seed and shared background input. Inflationary estimates are conditional on an effectively light minimally coupled spectator. Neither nonlinear domains nor anisotropies were evolved.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    results = Path(__file__).with_name('results')
    candidate = json.loads((results/'direct-subthreshold-friedmann.json').read_text())[0]
    match = json.loads((results/'direct-pole-matching-2026-11-25.json').read_text())[0]
    production = json.loads((results/'converged-independent-kernels-2026-11-25.json').read_text(encoding='utf-8-sig'))
    row = production['parent_angle240_sparse_k00685'][0]
    trigger = trigger_scan(candidate, row)
    args.output.write_text(json.dumps(dict(trigger=trigger), indent=2, allow_nan=False)+'\n')
    floquet = peak(.5)
    field = match['finite_field_screen']
    mediator = spectator_scan(np.sqrt(4*M*M-match['detuning_eV2']), field['lambda_S'], floquet)
    dark_scalar = spectator_scan(M, field['lambda_A_for_field_cap_point01_Mpl'], floquet)
    report = dict(trigger=trigger, mediator=mediator, light_scalar=dark_scalar,
                  numerical_note='Initial fixed amplitude bracket [1e13,1e22] eV missed the S roots. Replaced by a conformal-quartic-scaled bracket and checked the independent reference.')
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
