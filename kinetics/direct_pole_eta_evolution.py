"""Expanding metric-free linear eta/chi pilot on the saved background.

The two initial conditions are imposed at fixed z, one for chi and one
for an isotropic eta temperature perturbation. Optional early seeds use metric-free evolution
and an analytic angular monopole for ballistic transport. No adiabatic solution is imposed. Keep the collisionless eta
momentum and angle distribution instead of a fluid closure.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.special import roots_legendre

from direct_pole_eta_stability import MPC_CM
from direct_subthreshold_friedmann import solve_friedmann
from resonant_mediator_inverse import C, HBAR


def grids(order, angle_order):
    edges = np.array([0., .001, .01, .1, .4, 1., 2., 4., 8., 16., 32., 56.])
    u, w = roots_legendre(order)
    d = np.diff(edges)[:, None]/2
    q = ((edges[1:]+edges[:-1])[:, None]/2+d*u).ravel()
    dq = (d*w).ravel()
    u, w = roots_legendre(angle_order)
    occupations = 1/np.expm1(q)
    return q, dq, (u+1.)/2, w/2, occupations


def solve_mode(history, params, k, z_initial, order, angle_order, rtol, primordial_chi=False, primordial_eta=False):
    t = history['log_a']
    start = -np.log1p(z_initial)
    x_at = PchipInterpolator(t, history['x'])
    dx_at = PchipInterpolator(t, history['dx'])
    ln_h = PchipInterpolator(t, np.log(history['H']*HBAR))
    h_slope = ln_h.derivative()
    q, w, mu, wmu, f = grids(order, angle_order)
    dfdq = f*(1+f)
    vev = params['vev_eV']
    eta_mass0 = params['m_eta_at_vev_eV']
    g = eta_mass0/vev
    temp0 = params['Td_on_eV']/(params['z_on']+1)
    density_dual = params['Neta']/(2*np.pi**2)
    fm = (w*q*q*dfdq)[:, None]*wmu
    x0 = float(x_at(start))
    dx0 = float(dx_at(start))
    cw = params['eta_CW_over_X']
    mu2 = params['mu_eV']**2

    def local(t):
        x = float(x_at(t))
        h = np.exp(float(ln_h(t)))
        r = np.exp(-t)
        temp = temp0*r
        m_over_t = eta_mass0*x/temp
        e = np.hypot(q, m_over_t)
        beta = m_over_t/e
        m2_chi = mu2*(3*x*x-1 + 4*cw*(1+x*x*(3*np.log(max(1e-300, x*x))-1)))
        direct = density_dual*g*g*temp*temp*np.sum(w*q**4*f/e**3)
        wave = k*r*C*HBAR/(MPC_CM*h)
        return x, h, temp, e, beta, (m2_chi+direct)/h**2, wave

    early_start = -np.log(1e6)
    early_times = np.linspace(early_start, start, 2400) if primordial_eta else np.array([start])
    phase_times = np.r_[early_times, t[t > start]]
    physical_wave = k*np.exp(-phase_times-ln_h(phase_times))*C*HBAR/MPC_CM
    eta_m_over_t = eta_mass0*x_at(phase_times)*np.exp(phase_times)/temp0
    velocities = q[None, :]/np.hypot(q[None, :], eta_m_over_t[:, None])
    phase = cumulative_trapezoid(physical_wave[:, None]*velocities, phase_times,
                                 axis=0, initial=0.)
    ballistic = PchipInterpolator(phase_times, phase)

    def free_monopole(t):
        return q*np.sinc(ballistic(t)/np.pi)

    def free_feedback(t, beta, temp, h):
        return density_dual*g*temp**3/(vev*h*h)*(
            w*q*q*dfdq*beta @ free_monopole(t)
        )

    def rhs(t, packed):
        state = packed[:6].reshape(2, 3)
        response = packed[6:].reshape(q.size, angle_order, 2)
        u, v = state.real
        _, h, temp, e, beta, mass2, wave = local(t)
        susceptibility = fm*beta[:, None]*(density_dual*g*temp**3/(vev*h*h))
        feedback = np.sum(susceptibility[:, :, None]*response.real, axis=(0, 1))
        phidot = -(3+h_slope(t))*v - (wave*wave+mass2)*u
        phidot[:2] -= feedback
        phidot[1] -= free_feedback(t, beta, temp, h)
        streaming = (-1j*wave*mu[None, :])[:, :, None]
        delta_e_over_t = beta*eta_mass0/temp
        dr = streaming*(q[:, None, None]/e[:, None, None]*response
                        + delta_e_over_t[:, None, None]*u[None, None, :2])
        return np.concatenate((v, phidot, dr.ravel()))

    initial_chi = np.array([1., dx0/x0])
    chi_normalization = x0
    early_eta_mass_over_temperature = None
    initial_eta = np.zeros(2)
    if primordial_chi or primordial_eta:
        early_eval = np.linspace(early_start, start, 1200)
        early_eta_mass_over_temperature = float(np.max(
            abs(x_at(early_eval))*eta_mass0*np.exp(early_eval)/temp0
        ))

        def scalar_only(t, state):
            _, h, temp, _, beta, curvature, wave = local(t)
            y, dy = state.reshape(2, 2)
            acceleration = -(3+h_slope(t))*dy-(wave*wave+curvature)*y
            if primordial_eta:
                acceleration[1] -= free_feedback(t, beta, temp, h)
            return np.r_[dy, acceleration]

        early = solve_ivp(scalar_only, (early_start, start), [1., 0., 0., 0.],
                          method='DOP853', rtol=rtol, atol=rtol*.02, max_step=.018)
        if not early.success:
            raise RuntimeError(early.message)
        if primordial_chi:
            initial_chi = early.y[[0, 2], -1]
            chi_normalization = params['seed']
        if primordial_eta:
            initial_eta = early.y[[1, 3], -1]

    initial = np.zeros(6+2*q.size*angle_order, dtype=complex)
    initial[:6] = [initial_chi[0], initial_eta[0], initial_chi[0],
                   initial_chi[1], initial_eta[1], initial_chi[1]]
    tt = np.linspace(start, 0., 2200)
    sol = solve_ivp(rhs, (start, 0.), initial, t_eval=tt, method='DOP853',
                    rtol=rtol, atol=rtol*.02, max_step=.004)
    if not sol.success:
        raise RuntimeError(sol.message)
    field = x_at(tt)
    use = abs(field) >= .6
    physical_u = sol.y[:3].real*np.array([chi_normalization, 1., chi_normalization])[:, None]
    pmax = np.max(abs(physical_u[:, use]), axis=1)
    relative = np.max(abs(physical_u[:, use]/field[use]), axis=1)
    full = sol.y[6:].reshape(q.size, angle_order, 2, len(tt))
    number = (w*q*q*f).sum()
    eta_number = np.sum(fm[:, :, None, None]*full.real, axis=(0, 1))/number
    analytic_free = np.sinc(ballistic(tt)/np.pi) @ (w*q**3*dfdq)/number
    eta_number[1] += analytic_free
    threshold = np.flatnonzero(use & (abs(eta_number[1])*1e-5 >= 1.))
    z_nonlinear_eta = float(np.expm1(-tt[threshold[0]])) if threshold.size else None
    return dict(k_per_Mpc=k, z_initial=z_initial, radial_order=order, angle_order=angle_order,
                rtol=rtol, nfev=sol.nfev, trigger_at_initial=x0,
                primordial_chi=primordial_chi, primordial_eta=primordial_eta,
                ballistic_eta_number_contrast_at_coupling=float(analytic_free[0]),
                z_where_eta_number_reaches_one_for_1e_minus_5_input=z_nonlinear_eta,
                pre_coupling_max_eta_mass_over_temperature=early_eta_mass_over_temperature,
                max_abs_delta_x_after_roll=pmax.tolist(),
                max_abs_delta_x_relative_trigger_after_roll=relative.tolist(),
                today_delta_x=physical_u[:, -1].tolist(),
                max_abs_eta_number_fraction_after_roll=(np.max(abs(eta_number[:, use]), axis=1)
                                                       * np.array([chi_normalization, 1.])).tolist(),
                today_eta_number_fraction=(eta_number[:, -1]*np.array([chi_normalization, 1.])).tolist(),
                modes=['coupled initial fractional chi shift', 'coupled initial eta temperature shift',
                       'unforced chi control with identical chi shift'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', type=int, choices=(240, 320), required=True)
    parser.add_argument('--radial', type=int, default=7)
    parser.add_argument('--angle', type=int, default=12)
    parser.add_argument('--rtol', type=float, default=2e-7)
    parser.add_argument('--k', type=float, nargs='+', required=True)
    parser.add_argument('--z-initial', type=float, default=8.)
    parser.add_argument('--primordial-chi', action='store_true')
    parser.add_argument('--primordial-eta', action='store_true')
    parser.add_argument('--background-samples', type=int, default=22000)
    parser.add_argument('--seed-factor', type=float, default=1.)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source = Path(__file__).with_name('results')
    production = json.loads((source/'converged-independent-kernels-2026-11-25.json').read_text(encoding='utf-8-sig'))
    candidates = json.loads((source/'direct-subthreshold-friedmann.json').read_text())
    candidate = next(c for c in candidates if c['mesh'] == args.mesh)
    row = production[f'parent_angle{args.mesh}_sparse_k00685'][0]
    params = dict(candidate['params'])
    if args.seed_factor <= 0:
        parser.error('--seed-factor must be positive')
    params['seed'] *= args.seed_factor
    _, history = solve_friedmann(row, params, n=args.background_samples, rtol=1e-9,
                                 order=40, return_history=True)
    report = dict(production_mesh=args.mesh, background_samples=args.background_samples, seed_factor=args.seed_factor, modes=[],
                  scope='Local pilot on the prescribed expansion, no metric perturbations, no primordial adiabatic transfer, no nonlinear evolution. Values are per unit input amplitude in the three labeled modes.')
    for k in args.k:
        result = solve_mode(history, params, k, args.z_initial, args.radial, args.angle, args.rtol, args.primordial_chi, args.primordial_eta)
        report['modes'].append(result)
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
        print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
