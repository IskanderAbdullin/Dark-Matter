"""Archived-energy envelope for the specified isolated, late recycling network."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from archived_scalar_pair_reaction_reference import amplitude, points, snapshot_functions
from late_radiation_background_reference import TabulatedRadiation, MPL
from three_species_decay_onsager_reference import m, M, mp, y


def radiation_bound(eos, T, b_stop, safety=1.):
    R = TabulatedRadiation(eos, T)
    H_photon = np.pi / np.sqrt(45) * T*T/MPL
    lo = np.log(R.state(b_stop)[0])
    knots = R.rg.x
    assert np.array_equal(knots, R.sg.x)
    knots = np.r_[lo, knots[(knots > lo) & (knots < R.initial)], R.initial]
    lg, ls = R.rg(knots), R.sg(knots)
    # PCHIP does not overshoot between knots. Use independent worst-case
    # endpoints for rho and entropy; pad them outward by a tiny margin.
    sg0 = float(R.sg(R.initial))
    log_lower = .5*(np.minimum(lg[:-1], lg[1:]) - np.log(2)) + (
        2./3)*(sg0 - np.maximum(ls[:-1], ls[1:]))
    shell_bound = np.exp(log_lower - 1e-12)
    ta = (R.s0 - (3*knots + ls))/3
    assert np.all(np.diff(ta) < 0.) and abs(ta[-1]) < 1e-13
    assert abs(ta[0] - np.log(b_stop)) < 1e-12
    early = (np.diff(np.exp(-ta))/shell_bound).sum()
    late = 1/b_stop  # Assume a photon-only b^-2 lower floor after b_stop.
    independent = []
    for n in (5, 9):
        x, w = roots_legendre(n)
        d = np.diff(knots)[:, None]/2
        u = (knots[:-1, None]+knots[1:, None])/2 + d*x
        b = np.exp((R.s0-(3*u+R.sg(u)))/3)
        hh = np.exp(.5*(R.rg(u)-np.log(2)) + (2./3)*(sg0-R.sg(u)))
        jac = (3+R.sg.derivative()(u))/3
        assert np.min(jac) > 0.
        independent.append(float(np.sum((d*w)*jac/b/hh)))
    assert abs(independent[0]-independent[1]) < 1e-11
    assert independent[-1] < early
    return dict(floor_min=float(shell_bound.min()), photon_anchor_GeV=H_photon,
                b_stop=b_stop, panels=len(shell_bound), numerical_table_integral_times_Hph=independent,
                upper_integral_times_Hph=float((early+late)/safety),
                early_upper_times_Hph=float(early), late_photon_upper_times_Hph=late)


def moments(state):
    p = np.asarray(state["lab_momentum_edges_GeV"])
    f, _, g, _ = snapshot_functions(state)
    def sums(order):
        rows = []
        for mass, occ, deg in ((mp, f, 1), (m, g, 4)):
            _, E, u = points(p, mass, occ, order)
            rows.append([deg*u.sum(), deg*u@E])
        return np.asarray(rows)
    a, b = sums(32), sums(48)
    assert np.max(abs(a/b-1)) < 2e-10
    np.testing.assert_allclose(b[:, 0], np.asarray(state["represented_species_count_a3"]) /
                               state["scale_factor"]**3, rtol=2e-10, atol=0.)
    return b


def direct_spin_and_envelope_checks():
    z, I, one = np.zeros((2, 2)), np.eye(2), np.eye(4)
    pauli = (np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]),
             np.diag([1, -1]))
    gamma = [np.block([[I, z], [z, -I]])] + [
        np.block([[z, a], [-a, z]]) for a in pauli]

    def slash(v):
        return v[0]*gamma[0] - sum(c*g for c, g in zip(v[1:], gamma[1:]))

    err = []
    x, w = roots_legendre(96)
    max_ratio = max_geometric_ratio = 0.
    for s in np.geomspace(4*M*M*1.00001, 9., 30):
        E = np.sqrt(s)/2
        p = np.sqrt(s/4-M*M)
        k = np.sqrt(s/4-mp*mp)
        p1, p2 = [np.array([E, 0., 0., sign*p]) for sign in (1., -1.)]
        for a in np.linspace(-1., 1., 17):
            q = p1 - np.array([E, k*np.sqrt(1-a*a), 0., k*a])
            t = -q[1:] @ q[1:]
            X = slash(q) + m*one
            raw = np.trace((slash(p2)-M*one) @ X @ (slash(p1)+M*one) @ X).real/4
            scale = (m*m-t)*s
            analytical = scale/2-(M*M+t-mp*mp+2*M*m)**2/2
            err.append(abs(raw-analytical)/scale)
        # Final scalar beta cancels in this diagnostic ratio.
        ratio = (w @ amplitude(s, x, light=M, heavy=m, scalar=mp))/(4*s*y**4/(m*m))
        max_ratio = max(max_ratio, float(ratio))
        geom = np.sqrt(m**4 + m*m*(s-2*M*M-2*mp*mp) + (M*M-mp*mp)**2)
        max_geometric_ratio = max(max_geometric_ratio, float(ratio*geom/(m*m)))
    assert max(err)<2e-13 and max_ratio<1. and max_geometric_ratio<1.
    return dict(max_single_trace_defect_over_scale=float(max(err)),
                max_sampled_forward_cross_section_over_ceiling=max_ratio,
                max_sampled_forward_cross_section_over_geometric_ceiling=max_geometric_ratio)


def run(eos, paths, budget, goal=.120):
    K_inv = y**4/(np.pi*m*m)
    delta = M*M-m*m
    KD = y**4*(1+((M+m)**2-mp*mp)**2/delta**2)/(16*np.pi*M*M*np.sqrt(1-m*m/(M*M)))
    # Retain beta_phi in the forward ceiling; it cancels in detailed
    # balance. The remaining energy-weighted factor (4M^2/s)*beta_M
    # has its exact maximum at s=6M^2. K_inv already includes 8.
    K_phi_uniform = K_inv * 2/(3*np.sqrt(3))
    # The angle integral of (A-D cos(theta))^-1 is at most
    # 2/sqrt(A*A-D*D) by atanh(z) <= z/sqrt(1-z*z). The squared
    # denominator is increasing in s, with minimum at s=4*M*M.
    K_phi = K_phi_uniform * m*m/(m*m + M*M - mp*mp)
    K = max(K_phi, 2*KD)
    q = json.loads(Path(budget).read_bytes())["cases"]
    anchors = {int(r["case"].split("_")[1][1:]): r for r in q if r["case"].endswith("_12MeV")}
    rows = []
    ratio = mp/m
    for i, path in zip((16, 32), paths, strict=True):
        s = json.loads(Path(path).read_bytes())["snapshots"][-1]
        N, rho = moments(s).sum(axis=0)
        T = s["T_visible_GeV"]
        Y = anchors[i]["scalar_equivalent_inventory"]
        need = goal*(1+ratio)/2/Y
        out = []
        for cutoff,safety,energy_multiplier in ((100.,1.,1.),(500.,1.,1.),
                                                 (500.,.98,1.03),(500.,.98,1.05)):
            clk = radiation_bound(eos, T, cutoff, safety)
            gain = K*(rho*energy_multiplier)**2/(M*M*N)*clk["upper_integral_times_Hph"]/clk["photon_anchor_GeV"]
            coarse_gain = gain * max(K_phi_uniform, 2*KD)/K
            Ye = Y*(1+gain)
            phi = (Ye - ratio*goal)/(1-ratio)
            out.append(dict(clock=clk, H_multiplier=safety, energy_multiplier=energy_multiplier,
                            upper_added_count_per_initial_represented=float(gain),
                            uniform_angle_baseline_upper_inventory_multiplier=float(1+coarse_gain),
                            upper_inventory_multiplier=float(1+gain),
                            upper_scalar_Omega_at_fixed_total=float(phi),
                            half_at_reference_excluded=bool(1+gain < need)))
        rows.append(dict(q=i, T_GeV=T, represented_number_GeV3=float(N),
                         represented_energy_GeV4=float(rho), conservative_initial_Y=Y,
                         multiplier_required_for_half=need, bounds=out))
    return dict(scope="Isolated homogeneous represented system, no future external source "
                      "or reentry, stated vacuum channels only. Allows arbitrary elastic "
                      "redistribution and Bose stimulation inside those channels; "
                      "conditioned visible clock, positive pressures.",
                target_total=goal,K_phi_GeVminus2=float(K_phi),
                K_phi_uniform_GeVminus2=float(K_phi_uniform),
                combined_K_GeVminus2=float(K),
                twice_stable_KD_GeVminus2=float(2*KD),
                independent_spin_checks=direct_spin_and_envelope_checks(),cases=rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("eos")
    parser.add_argument("budget")
    parser.add_argument("paths", nargs=2)
    a = parser.parse_args()
    print(json.dumps(run(a.eos, a.paths, a.budget), indent=2, allow_nan=False))
