"""Local collisionless eta response to the mass-setting trigger field.

Freeze a homogeneous background snapshot, keep relativistic velocities and
the decoupled Bose momentum distribution, and solve the growing scalar /
Vlasov dispersion relation. Includes scalar Hubble friction as a local
estimate. This is not an expanding perturbation integration.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid
from scipy.optimize import brentq
from scipy.special import roots_legendre

from direct_subthreshold_friedmann import solve_friedmann
from resonant_mediator_inverse import C, HBAR

MPC_CM = 3.085677581e24


def angular_response(y):
    """1 - atan(y)/y, where y = k*v/gamma."""
    y = np.asarray(y)
    low = y < .003
    ans = np.empty_like(y)
    u = y[low]**2
    ans[low] = u*(1/3+u*(-1/5+u/7))
    ans[~low] = 1-np.arctan(y[~low])/y[~low]
    return ans


class EtaResponse:
    def __init__(self, order, qmin):
        edges = np.r_[0., np.geomspace(qmin, .1, 15), .3, 1., 2., 4., 8., 16., 32., 56.]
        x, w = roots_legendre(order)
        radius = np.diff(edges)/2
        self.q = ((edges[1:]+edges[:-1])[:, None]/2+radius[:, None]*x).ravel()
        self.weight = (radius[:, None]*w).ravel()
        self.f = 1/np.expm1(self.q)
        self.minus_df = self.f*(1+self.f)

    def coefficients(self, mass, temperature, coupling, species):
        energy = np.hypot(self.q, mass/temperature)
        velocity = self.q/energy
        g = species*coupling**2/(2*np.pi**2)
        direct = g*temperature**2 * np.sum(self.weight*self.q**4/energy**3*self.f)
        susceptibility = g*mass**2 * self.weight*velocity*self.minus_df
        return direct, velocity, susceptibility

    def rate(self, k, H, effective_mass2, velocity, susceptibility):
        """Return gamma/H above the finite-growth reporting floor."""
        k_h = k/H
        k_v = k_h*velocity
        curvature = effective_mass2/H**2
        s = susceptibility/H**2

        def dispersion(g):
            return g*(g+3) + k_h**2 + curvature - s @ angular_response(k_v/g)

        floor = 1e-5
        if dispersion(floor) >= 0:
            return 0.
        high = max(1., np.sqrt(max(0., s.sum()-k_h**2-curvature))+1)
        return brentq(dispersion, floor, high, xtol=1e-8, rtol=1e-9)


def solve(row, candidate, order, qmin, snapshots, nk):
    params = candidate['params']
    _, history = solve_friedmann(row, params, n=11000, rtol=2e-9,
                                 order=40, return_history=True)
    time = history['log_a']
    criteria = (np.abs(history['x'])**4 >= .5)
    first = np.flatnonzero(criteria)[0]
    t = np.linspace(time[first], 0., snapshots)
    r = np.exp(-t)
    field, dx, hubble = (np.interp(t, time, history[k]) for k in ('x', 'dx', 'H'))
    hubble *= HBAR
    temperature = params['Td_on_eV']/(params['z_on']+1)*r
    coupling = params['m_eta_at_vev_eV']/params['vev_eV']
    mass = params['m_eta_at_vev_eV']*np.abs(field)
    k = np.geomspace(.001, 30., nk)
    grid = np.zeros((snapshots, nk))
    m2 = np.zeros_like(t)
    cold_scale = np.zeros_like(t)
    eta = EtaResponse(order, qmin)
    for i, (x, tt, mm, hh, rr) in enumerate(zip(field, temperature, mass, hubble, r)):
        cw = params['eta_CW_over_X']
        potential = params['mu_eV']**2*(3*x*x-1 + 4*cw*(1 + x*x*(3*np.log(x*x)-1)))
        direct, velocities, susceptibility = eta.coefficients(mm, tt, coupling, params['Neta'])
        m2[i] = potential+direct
        cold_scale[i] = _cold_scale(coupling, params['Neta'], tt, mm)/hh
        for j, comoving in enumerate(k):
            grid[i, j] = eta.rate(comoving*rr*C*HBAR/MPC_CM, hh, m2[i], velocities, susceptibility)
    k_h = k[None, :]*r[:, None]*C*HBAR/MPC_CM/hubble[:, None]
    # These necessary scale-separation cuts do not guarantee local validity.
    usable = (grid >= 5.) & (grid >= 3*abs(dx/field)[:, None]) & (k_h >= 10.)
    unforced = np.maximum((np.sqrt(np.maximum(9-4*(k_h**2+(m2/hubble**2)[:, None]), 0.))-3)/2, 0.)
    integrated_unforced = trapezoid(unforced, t, axis=0)
    integrated = trapezoid(grid, t, axis=0)
    integrated_fast = trapezoid(np.where(usable, grid, 0.), t, axis=0)
    i, j = np.unravel_index(np.argmax(grid), grid.shape)
    j_best = int(np.argmax(integrated_fast))
    return dict(mesh=candidate['mesh'], momentum_quad_order=order, smallest_q_edge=qmin,
                max_local_gamma_over_H=float(grid[i, j]), z_at_max_local=float(r[i]-1),
                k_at_max_local_per_Mpc=float(k[j]), z_start=float(r[0]-1),
                max_integrated_frozen_gamma=float(max(integrated)),
                max_integrated_fast_subset=float(max(integrated_fast)),
                k_at_max_integrated_fast_per_Mpc=float(k[j_best]),
                minimum_effective_scalar_m2_over_H2=float(min(m2/hubble**2)),
                maximum_cold_limit_gamma_scale_over_H=float(max(cold_scale)),
                max_integrated_unforced_scalar_gamma=float(max(integrated_unforced)),
                maximum_fixed_snapshot_eta_increment=float(max(integrated-integrated_unforced)),
                k_per_Mpc=k.tolist(), z=(r-1).tolist(), trigger=field.tolist(),
                effective_scalar_m2_over_H2=(m2/hubble**2).tolist(),
                local_peak_gamma_over_H_by_z=np.max(grid, axis=1).tolist(),
                integrated_gamma_by_k=integrated.tolist(), fast_integrated_gamma_by_k=integrated_fast.tolist(),
                gamma_at_max_fast_k_over_H=grid[:, j_best].tolist(),
                mass_eV=mass.tolist(), temperature_eV=temperature.tolist(),
                scope='Frozen homogeneous snapshots, relativistic collisionless eta response, no metric perturbations, expanding Liouville history, or nonlinear eta evolution.')


def _cold_scale(coupling, species, temperature, mass):
    # sqrt(g^2*n/m). zeta(3) sets the frozen massless thermal number.
    return np.sqrt(coupling**2*species*1.202056903159594*temperature**3/(np.pi**2*mass))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).with_name('results')
    candidates = json.loads((root/'direct-subthreshold-friedmann.json').read_text())
    source = json.loads((root/'converged-independent-kernels-2026-11-25.json').read_text(encoding='utf-8-sig'))
    reports = []
    for candidate, quad, qmin, ns, nk in zip(candidates, (16, 28), (1e-11, 1e-14), (336, 672), (104, 168)):
        row = source[f'parent_angle{candidate["mesh"]}_sparse_k00685'][0]
        result = solve(row, candidate, quad, qmin, ns, nk)
        reports.append(result)
        args.output.write_text(json.dumps(reports, indent=2, allow_nan=False)+'\n')
        print(json.dumps({k: v for k, v in result.items() if not isinstance(v, list)}), flush=True)


if __name__ == '__main__':
    main()
