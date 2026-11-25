"""Time-ordered homogeneous NR Bose scattering through a virtual s-channel pole.

The exact NR angle average multiplies the existing conservative nodal graph.
The Friedmann history determines the physical velocity of a comoving pair.
Only AA <-> AA is evolved. Crossed terms, relativistic corrections, medium
corrections, other species and inhomogeneities are outside this evolution.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import BDF, cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.special import roots_legendre

from direct_subthreshold_activation import VREF
from direct_subthreshold_friedmann import solve_friedmann
from elastic_implicit_trial import ImplicitElastic, merge_grid
from resonant_mediator_inverse import C
from velocity_only_exposure import RHO0


def pole_average(total_p2, low, high, detuning):
    """Average [b/(b+2 sum(p_i^2)-P^2)]^2 over the allowed interval dP."""
    d = detuning + total_p2
    root = np.sqrt(d)
    # Rational primitive difference factored by the interval width. The
    # logarithmic difference uses log1p to stay accurate in the contact limit.
    width = high - low
    den_high = detuning + np.maximum(total_p2 - high**2, 0.)
    den_low = detuning + np.maximum(total_p2 - low**2, 0.)
    log_part = np.log1p(width / (root + low)) + np.log1p(
        width * (root + high) / den_high
    )
    return detuning**2 / (2*d) * (
        (d + high*low) / (den_high*den_low) + log_part / (2*root*width)
    )


def event_geometry(model):
    p1, p2, p3 = model.p[model.idx[:3]]
    p4 = np.sqrt(np.maximum(p1**2 + p2**2 - p3**2, 0.))
    low = np.maximum(abs(p1-p2), abs(p3-p4))
    high = np.minimum(p1+p2, p3+p4)
    # Use the same representative four-momentum geometry as the base graph.
    assert np.all(high > low)
    return 2*(p1**2 + p2**2), low, high


def exposure(candidate, row):
    settings = candidate['full_friedmann']
    summary, history = solve_friedmann(
        row, candidate['params'], rtol=settings['rtol'], n=settings['n'],
        order=settings['order'], return_history=True
    )
    late = history['r'] <= 101*(1+1e-14)
    a, r, field, hubble = (history[k][late] for k in ('log_a', 'r', 'x', 'H'))
    integrated = cumulative_trapezoid(r**4*field**4/hubble, a, initial=0.)
    tau = candidate['sigma0']*RHO0*summary['v0']*integrated[-1]
    unique = np.r_[True, np.diff(integrated) > 0]
    return float(tau), summary['v0'], PchipInterpolator(
        integrated[unique]/integrated[-1], np.log(r[unique])
    )


def averaged_weight(geometry, log_redshift, today_detuning, start, stop, nodes):
    points, weights = roots_legendre(nodes)
    redshifts = log_redshift((stop + start + (stop-start)*points)/2)
    answer = np.zeros_like(geometry[0])
    for weight, log_r in zip(weights, redshifts):
        answer += weight/2 * pole_average(*geometry, today_detuning*np.exp(-2*log_r))
    return answer


def run(row, candidate, n, slices, nodes, rtol=1e-6, contact=False):
    f0 = row['frozen_cell_spectrum']
    peak = row['final']['peak_bin_q']
    edges, f = merge_grid(f0['edges'], f0['f_phi'], n, qmin=0., qmax=f0['edges'][-1])
    edges = edges/peak
    # Stable limit: dK_rel*(2*m) -> d(p^2), including all three energy measures.
    model = ImplicitElastic(edges, 1e10)
    base = model.rate*(2*model.m)**3
    geometry = event_geometry(model)
    tau, v0, log_r = exposure(candidate, row)
    detuning = candidate['a']*(VREF/v0)**2
    vol = model.vol
    number = vol @ f
    energy = (vol*model.K) @ f
    energy_measure = vol*model.p**2
    moments = [vol*model.p**4]
    for cut in (2., 6., 20.):
        lo, hi = edges[:-1], edges[1:]
        moments.append(vol*np.maximum(hi**3 - np.maximum(lo, cut)**3, 0.)/(hi**3-lo**3))
    moments = np.asarray(moments)
    reference = moments @ f
    norm = np.maximum(f, 1.)
    initial = np.maximum(f, 1e-38)
    state = initial/norm
    scale = tau/(128*np.pi**4*number)
    max_negative = 0.
    steps = []
    cuts = [0.01, 0.03, 0.1]
    occupancy_speed = edges[1:]*v0*101/C
    early_tail = {str(c): [float(vol[occupancy_speed>c]@f[occupancy_speed>c]/number),
                           float(energy_measure[occupancy_speed>c]@f[occupancy_speed>c]/(energy_measure@f))]
                  for c in cuts}
    eq = 1/np.expm1(.3 + 80*model.K/model.K[-1])
    equilibrium_residual = 0.
    time_edges = np.linspace(0., 1., slices+1)
    for lo, hi in zip(time_edges[:-1], time_edges[1:]):
        weight = 1. if contact else averaged_weight(geometry, log_r, detuning, lo, hi, nodes)
        model.rate = base * weight

        def rhs(_, y):
            return model.rates_jac(y*norm, False)*scale/norm

        def jac(_, y):
            return model.rates_jac(y*norm, True)*scale*norm[None, :]/norm[:, None]

        if not steps:
            eq_test = model.rates_jac(eq, False)
            flux_scale = np.sum(model.rate*eq[model.i]*eq[model.j])
            equilibrium_residual = float(vol @ abs(eq_test)/flux_scale)
        sol = BDF(rhs, lo, state, hi, jac=jac, rtol=rtol, atol=1e-13, max_step=(hi-lo)/8)
        while sol.status == 'running':
            note = sol.step()
            max_negative = max(max_negative, vol @ np.maximum(-sol.y*norm, 0.)/number)
            if sol.status == 'failed' or max_negative > 1e-9:
                raise RuntimeError(f'BDF stopped at exposure fraction {sol.t}: {note}; {max_negative=}')
        state = sol.y
        g = state*norm
        steps.append(dict(log_a_start=float(-log_r(lo)), log_a_stop=float(-log_r(hi)),
                          nfev=sol.nfev, moment_changes=(moments @ g/reference-1).tolist()))
    return dict(n=model.n, events=len(base), slices=slices, nodes=nodes, rtol=rtol, contact=contact,
                tau=tau, v0_cm_s=v0, today_detuning=detuning,
                maximum_grid_speed_over_c_z100=float(occupancy_speed[-1]),
                initial_nonrelativistic_diagnostics_z100=early_tail,
                equilibrium_residual_over_dilute_flux=equilibrium_residual,
                number_residual=float(vol @ g/number-1), energy_residual=float((vol*model.K) @ g/energy-1),
                maximum_negative_number_fraction=float(max_negative),
                initial_floor_number_fraction=float(vol @ (initial-f)/number),
                moments=['q4', 'number_above_2_peak', 'number_above_6_peak', 'number_above_20_peak'],
                final_moment_changes=(moments @ g/reference-1).tolist(), steps=steps,
                edges_over_peak=edges.tolist(), initial_f=f.tolist(), final_f=g.tolist())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', type=int, choices=(240, 320), required=True)
    parser.add_argument('--n', type=int, default=128)
    parser.add_argument('--slices', type=int, default=4)
    parser.add_argument('--nodes', type=int, default=16)
    parser.add_argument('--rtol', type=float, default=1e-6)
    parser.add_argument('--contact', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    results = Path(__file__).with_name('results')
    source = json.loads((results/'converged-independent-kernels-2026-11-25.json').read_text(encoding='utf-8-sig'))
    candidates = json.loads((results/'direct-subthreshold-friedmann.json').read_text())
    candidate = next(c for c in candidates if c['mesh'] == args.mesh)
    row = source[f'parent_angle{args.mesh}_sparse_k00685'][0]
    out = run(row, candidate, args.n, args.slices, args.nodes, args.rtol, args.contact)
    out['production_mesh'] = args.mesh
    args.output.write_text(json.dumps(out, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k: v for k, v in out.items() if k not in ('initial_f', 'final_f', 'edges_over_peak', 'steps')}))


if __name__ == '__main__':
    main()
