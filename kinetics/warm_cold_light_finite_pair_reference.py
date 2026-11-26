"""Fixed tagged-momentum scalar pair collision, averaged over one native cell."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import expit, roots_legendre

from archived_scalar_pair_reaction_reference import amplitude, gross
from frw_exterior_restart_reference import import_restart
from joint_visible_heavy_adiabatic_reference import snapshot_functions
from joint_prompt_two_moment_evolution_reference import Cells
from late_radiation_background_reference import MPL, TabulatedRadiation
from three_species_decay_onsager_reference import m, mp
from warm_cold_light_scalar_pair_reference import calc as at_rest
from warm_light_occupation_stage_audit_reference import digest, read
from warm_prompt_reentry_shape_audit_reference import get_mom


def tagged(pt, F, B, fe, be, order, batch=24):
    """Return gain and loss; bath occupations outside represented support are zero.

    Tag a single spin and charge; partner has two spin states. Initial direction
    is isotropic. The identical scalar factor is in the invariant prefactor.
    Both final scalars cannot be occupied beyond the sum of their support limits.
    Vacuum final states are still available for the forward process.
    """
    nr, ni, no, na = order
    Et = np.hypot(m, pt)
    Emax = max(fe[-1], 2 * be[-1] - Et)
    edge = np.unique(np.clip(np.r_[m, fe, 2 * be - Et, Emax], m, Emax))
    x, w = roots_legendre(nr)
    u, wi = roots_legendre(ni)
    z, wo = roots_legendre(no)
    az = np.cos(np.pi * (np.arange(na) + .5) / na)
    h = np.diff(edge) / 2
    es = (((edge[:-1] + edge[1:]) / 2)[:, None] + h[:, None] * x).ravel()
    ew = (h[:, None] * w).ravel()
    f0 = F(np.array([Et]))[0]
    totals = np.zeros(2)
    worst = 0.
    for j in range(0, len(es), batch):
        E = es[j:j + batch, None]
        p = np.sqrt((E - m) * (E + m))
        s = 2 * m * m + 2 * (Et * E - pt * p * u)
        ci = np.sqrt(np.maximum(s / 4 - m * m, 0.))
        cf = np.sqrt(s / 4 - mp * mp)
        boost = np.sqrt((p - pt)**2 + 2 * p * pt * (1 + u)) / np.sqrt(s)
        D = boost * cf
        mid = np.broadcast_to((Et + E) / 2, s.shape)
        # Swapping the outgoing scalars changes neither occupation product
        # nor the azimuth-averaged amplitude. Use the positive cosine only.
        cut = np.clip(abs((be - mid[..., None]) / np.maximum(D[..., None], 1e-280)), 0., 1.)
        cut = np.sort(np.concatenate([np.zeros_like(cut[..., :1]), cut,
                                      np.ones_like(cut[..., :1])], axis=-1), axis=-1)
        d = np.diff(cut, axis=-1) / 2
        i, k, l = np.where(d > 0.)
        co = (cut[i, k, l, None] + cut[i, k, l + 1, None]) / 2 + d[i, k, l, None] * z
        aw = 2 * d[i, k, l, None] * wo
        c0 = np.clip((Et - E) / (2 * np.maximum(boost * ci, 1e-280)), -1., 1.)[i, k, None]
        a0 = c0 * co
        a1 = np.sqrt(np.maximum(1 - c0 * c0, 0.)) * np.sqrt(1 - co * co)
        a = np.clip(a0[..., None] + a1[..., None] * az, -1., 1.)
        amp = amplitude(s[i, k, None, None], a).mean(axis=-1)
        f1 = F(E[:, 0])[i, None]
        b0 = B(mid[i, k, None] + D[i, k, None] * co)
        b1 = B(mid[i, k, None] - D[i, k, None] * co)
        G = (1 - f0) * (1 - f1) * b0 * b1
        L = f0 * f1 * (1 + b0) * (1 + b1)
        prior = p * ew[j:j + batch, None] / np.pi**2 * (wi / 2) * cf / (64 * np.pi * Et * np.sqrt(s))
        Q = prior[i, k, None] * aw * amp
        totals += [np.sum(Q * G), np.sum(Q * L)]
        worst = max(worst, np.max(abs(G - L) / np.maximum(np.maximum(G, L), 1e-280)))
    return totals, float(worst)


def controls():
    T, mu = .022, -.02
    B = lambda E: 1 / np.expm1((E - mu) / T)
    F = lambda E: expit((mu - E) / T)
    be = np.hypot(np.linspace(0., 1., 15), mp)
    fe = np.hypot(np.linspace(0., 1., 15), m)
    thermal = []
    for p in (0., .014, .085):
        v, d = tagged(p, F, B, fe, be, (10, 6, 10, 8))
        defect = float(abs(v[0] / v[1] - 1))
        assert d < 6e-14 and defect < 1e-14
        thermal.append(dict(p_GeV=p, event_balance=d, integrated_balance=defect))
    return thermal


def one(eos, S, orders):
    p = np.asarray(S["lab_momentum_edges_GeV"])
    B, be, F, fe = snapshot_functions(S)
    rho = TabulatedRadiation(eos, S["T_visible_GeV"]).state(1.)[1]
    grav = get_mom(p, [S["entropy_variable_scalar"], S["entropy_variable_fermion"]])[2].sum()
    H = np.sqrt((rho + grav + import_restart(S["positive_exterior_restart"]).energy(1.)) / 3) / MPL
    ref, _ = at_rest(B, be, F, fe, 28, 24)
    result = []
    for a, r, i, o, z in orders:
        t, w = roots_legendre(a)
        x = p[1] * (1 + t) / 2
        wx = w * x**2
        fixed = np.array([tagged(v, F, B, fe, be, (r, i, o, z))[0] for v in x])
        vals = wx @ fixed / (wx.sum() * H)
        at0, _ = tagged(0., F, B, fe, be, (r, i, o, z))
        defect = float(np.max(abs(at0 / ref - 1.)))
        result.append(dict(order=[a, r, i, o, z], first_cell_gain_loss_per_H=vals.tolist(),
                           first_cell_net_per_H=float(vals[0] - vals[1]), rest_match=defect))
    assert result[-1]["rest_match"] < 3e-7
    assert abs(result[-1]["first_cell_net_per_H"] / result[-2]["first_cell_net_per_H"] - 1) < 2e-5
    native = []
    for order in ((3, 7, 6, 2), (4, 10, 10, 4)):
        _, _, gain = gross(p, mp, m, B, F, True, *order, project=True)
        _, _, loss = gross(p, m, mp, F, B, False, *order, project=True)
        norm = H * Cells(p, m, 4).vol[0]
        native.append(dict(order=order, first_cell_gain_loss_per_H=[float(gain[1, 0, 1] / norm),
                                                                     float(loss[0, 0, 1] / norm)]))
    v = np.array([z["first_cell_gain_loss_per_H"] for z in native])
    direct = np.array(result[-1]["first_cell_gain_loss_per_H"])
    assert np.max(abs(v[:, 1] / direct[1] - 1.)) < 2e-6
    assert v[0, 0] < 2e-5 * direct[0] and .5 < v[1, 0] / direct[0] < .55
    return dict(T_GeV=S["T_visible_GeV"], first_edge_GeV=float(p[1]), rest_gain_loss_per_H=(ref / H).tolist(),
                native_gross=native, rows=result)


def run(eos, D):
    budget = read(D / "warm-recoil-endpoint-channel-budget-2026-11-26.json")
    cases = []
    for q, native in zip((32, 64), budget["cases"], strict=True):
        entry = read(D / f"warm-recoil-crossing-q{q}-terminal-2026-11-26.json.gz")[-1]
        S = entry["snapshot"]
        assert digest(S) == entry["snapshot_sha256"]
        a = one(eos, S, ((4, 8, 6, 8, 8), (6, 12, 9, 12, 12), (8, 18, 12, 18, 16)))
        point = a["rows"][-1]["first_cell_net_per_H"]
        b = [float(z["first_cold_shell_rates_per_log_a"][2] / point) for z in native["controls"]]
        for frozen, z in zip(a["native_gross"], native["controls"], strict=True):
            gain, loss = frozen["first_cell_gain_loss_per_H"]
            assert abs((gain - loss) / z["first_cold_shell_rates_per_log_a"][2] - 1.) < 2e-12
        cases.append(dict(q=q, sha256=digest(S), native_to_independent_net_ratio=b, **a))
    return dict(scope="Frozen retained bath on two actual post-recoil checkpoints; no correction to evolved states.",
                thermal=controls(), cases=cases)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos", type=Path)
    p.add_argument("--tapes", type=Path, default=Path("results"))
    a = p.parse_args()
    print(json.dumps(run(a.eos, a.tapes), indent=2, allow_nan=False))
