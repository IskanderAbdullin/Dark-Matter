"""Compare cut counts with excess energy on the same saved prompt phases."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from joint_visible_heavy_adiabatic_reference import angle_panels, snapshot_functions
from three_species_decay_onsager_reference import Es, M, P, m, mp, width
from warm_light_occupation_stage_audit_reference import digest, read
from warm_prompt_reentry_shape_audit_reference import get_mom


def resident(S, cut, order):
    fs, _, ff, _ = snapshot_functions(S)
    edges = np.asarray(S["lab_momentum_edges_GeV"])
    x, w = roots_legendre(order)
    out = []
    for mass, g, f in ((mp, 1, fs), (m, 4, ff)):
        qcut = np.sqrt((cut - mass) * (cut + mass))
        p = np.unique(np.r_[qcut, edges[edges > qcut]])
        half = np.diff(p)[:, None] / 2
        q = p[:-1, None] + half * (1 + x)
        energy = np.hypot(q, mass)
        dN = g / (2 * np.pi**2) * w * half * q**2 * f(energy)
        out.append([float(dN.sum()), float((dN * (energy - cut)).sum()),
                    float((dN * energy).sum())])
    return np.asarray(out)


def currents(S, cut, order):
    fs, es, ff, ef = snapshot_functions(S)
    A, B = np.unique(np.r_[es, cut]), np.unique(np.r_[ef, cut])

    def parent_walls(walls, mass, rest):
        q = np.sqrt((walls - mass) * (walls + mass))
        hi = walls * rest + q * P
        return np.r_[M * (rest**2 + q**2) / hi, M * hi / mass**2]

    lo = max(cut, M)
    hi = es[-1] + ef[-1]
    boundaries = np.unique(np.clip(np.r_[lo, hi, parent_walls(A, mp, Es),
        parent_walls(B, m, M-Es), (A[:, None] + B).ravel()], lo, hi))
    x, w = roots_legendre(order)
    half = np.diff(boundaries)[:, None] / 2
    E = (boundaries[:-1, None] + half * (1 + x)).ravel()
    dE = (half * w).ravel()
    out = np.zeros((2, 2))
    abs_ledgers = np.zeros_like(out)
    for i in range(0, len(E), 96):
        ep = E[i:i+96]
        s, W = angle_panels(ep, A, B, angle=order)
        f = ep[:, None, None] - s
        W *= ((s >= es[0]) & (s <= es[-1]) & (f >= ef[0]) & (f <= ef[-1]))
        b, c = fs(s), ff(f)
        inv = W * b * c
        loss = W * (1+b) * (1-c)
        formation = inv.sum(axis=(1, 2))
        norm = formation + loss.sum(axis=(1, 2))
        h = np.divide(formation, norm, where=norm>0, out=np.zeros_like(norm))[:, None, None]
        release, consume = h * loss, (1-h) * inv
        measure = dE[i:i+96] * 2 / np.pi**2 * width * M * np.sqrt((ep-M)*(ep+M))
        for j, energy in enumerate((s, f)):
            for k, weight in enumerate((energy > cut, np.maximum(energy-cut, 0.))):
                out[j, k] += measure @ ((release-consume)*weight).sum(axis=(1, 2))
                abs_ledgers[j, k] += measure @ ((release+consume)*weight).sum(axis=(1, 2))
    return out, abs_ledgers


def run(folders, steps, orders=(6, 10)):
    out = []
    for folder in folders:
        plan = read(folder/"summary.json")
        H = plan["H_GeV"]
        for step in steps:
            a = read(folder/f"step-{step}.json")
            S = a["snapshot"]
            assert digest(S) == a["hash"] and a["plan"] == plan["plan"]
            q = np.asarray(S["lab_momentum_edges_GeV"])
            Ns, _, Eshell = get_mom(q, [S["entropy_variable_"+s] for s in ("scalar", "fermion")])
            N = Ns.sum(axis=1)
            rho = Eshell.sum()
            cut = a["plan"]["normalized_cut_GeV"]
            R = resident(S, cut, 16)
            np.testing.assert_allclose(R[:, 2], cut*R[:, 0] + R[:, 1], rtol=5e-15)
            np.testing.assert_allclose(resident(S, cut, 24), R, rtol=2e-11)
            norm = np.stack((N, np.full(2, rho)), axis=1)
            rates = []
            for order in orders:
                signed, gross = currents(S, cut, order)
                signed, gross = signed/(norm*H), gross/(norm*H)
                assert abs(signed[:, 1].sum()-a["signed_native_flux"]) < 6e-5*gross[:, 1].sum()
                assert np.all(abs(signed) <= gross*(1+3e-15))
                rates.append(dict(order=order, signed=signed.tolist(), gross=gross.tolist()))
            diff = abs(np.asarray(rates[-1]["signed"])-rates[0]["signed"])
            assert np.all(diff < 6e-5*np.asarray(rates[-1]["gross"]))
            out.append(dict(hard_GeV=a["plan"]["hard_GeV"],halfwidth_GeV=a["plan"]["shell_halfwidth_GeV"],
                Hdt=a["time"], sha256=a["hash"], N_above_over_species_N=(R[:, 0]/N).tolist(),
                W_over_rho=(R[:, 1]/rho).tolist(), full_energy_above_over_rho=(R[:, 2]/rho).tolist(),
                current_controls=rates))
    return dict(scope="Independent angular count and convex-energy rates at saved fixed-background "
                      "phases; no counting of parent or other processes.",
                species=["scalar","light"], columns=["N_above / full species N / H-time",
                   "convex excess energy / full dark energy / H-time"], cases=out)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("folders", type=Path, nargs="+")
    p.add_argument("--steps", type=int, nargs="+", default=[0, 5, 6])
    args = p.parse_args()
    print(json.dumps(run(args.folders, args.steps), indent=2, allow_nan=False))
