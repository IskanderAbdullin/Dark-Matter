"""Hash, ledger, and shared-support audit of own-history continuation after recoil."""

import argparse
import json
from pathlib import Path

import numpy as np

from frw_exterior_restart_reference import import_restart
from joint_prompt_two_moment_evolution_reference import Collision
from late_radiation_background_reference import MPL, TabulatedRadiation
from three_species_decay_onsager_reference import m, mp
from warm_coupled_serial_restart_reference import sha
from warm_history_lower_matched_audit_reference import measures
from warm_light_occupation_stage_audit_reference import digest, read
from warm_near_rest_light_inverse_decay_reference import window


def check(eos, seed, folder, cap):
    t = read(seed)
    assert all(digest(row["snapshot"]) == row["snapshot_sha256"] for row in t)
    S = t[-1]["snapshot"]
    A = S["scale_factor"]
    E0 = measures(S)
    M0 = E0[0].sum() * A**3
    ext0 = import_restart(S["positive_exterior_restart"]).array()
    base_count = ext0.shape[1]
    last = digest(S)
    dh = sha(seed)
    eh = sha(eos)
    plan = None
    V = [E0 * A**3]
    rows = []
    for i in range(cap + 1):
        if i:
            C = read(folder / f"step-{i}.json")
            assert C["index"] == i and C["previous"] == last
            assert C["plan"]["tape_sha256"] == dh and C["plan"]["EOS_sha256"] == eh
            plan = C["plan"] if plan is None else plan
            assert plan == C["plan"]
            S = C["snapshot"]
            last = digest(S)
            assert C["snapshot_sha256"] == last
        a = S["scale_factor"]
        b = a / A
        assert abs(b / (1.10**(i / 24)) - 1.) < 1e-13
        p = np.array(S["lab_momentum_edges_GeV"])
        Q = measures(S)
        ext = import_restart(S["positive_exterior_restart"])
        ar = ext.array()
        restored = ar[:, :base_count].copy()
        restored[[0, 5]] *= np.array([b, b**3])[:, None]
        np.testing.assert_allclose(restored, ext0, rtol=1.5e-13, atol=0.)
        num = float(abs((Q[0].sum() + ar[5, base_count:].sum()) * a**3 / M0 - 1.))
        assert num < 6e-13
        T = S["T_visible_GeV"]
        rho = TabulatedRadiation(eos, T).state(1.)[1]
        H = np.sqrt((rho + Q[1].sum() + ext.energy(1.)) / 3) / MPL
        c = Collision(p)
        th = [np.array(S["entropy_variable_" + j]) for j in ("scalar", "fermion")]
        native = float(c.rate(th)[1, 0, 0] / (H * c.c[1].vol[0]))
        ecap = np.hypot(mp, p[-1])
        if ecap < window(p[1])[-2]:
            assert native == 0. and not np.any(c.J == 0)
        f = Q[0, 1, 0] / c.c[1].vol[0]
        slope = None if not i else (f - rows[-1]["cold_mean_occupation"]) / np.log(b / rows[-1]["relative_a"])
        rows.append(dict(step=i, sha256=last, T_GeV=T, relative_a=b, cumulative_number_residual=num,
                         local_energy_defect=C["operator_summary"]["energy_ledger"] if i else None,
                         newton_history_only=C["audit"]["solver"]["used_history"] if i else None,
                         scalar_edge_GeV=float(ecap), native_events=int(np.sum(c.J == 0)),
                         selected_prompt_cold_net_per_H=native, cold_mean_occupation=float(f),
                         accepted_cold_slope=None if slope is None else float(slope)))
        if i:
            V.append(Q * a**3)
    return rows, np.array(V), p


def run(eos, work, n):
    r, x, p = check(eos, work / "seed32.json", work / "q32", n)
    s, y, P = check(eos, work / "seed64.json", work / "q64", n)
    np.testing.assert_allclose(p, P[:len(p)], rtol=1e-13, atol=0.)
    np.testing.assert_allclose([j["T_GeV"] for j in r], [j["T_GeV"] for j in s], rtol=1e-13, atol=0.)
    assert all(j["selected_prompt_cold_net_per_H"] == 0 for j in r)
    assert all(j["selected_prompt_cold_net_per_H"] < 0 for j in s)
    den = np.vstack((x[0, 0].sum(axis=1), np.full((2, 2), x[0, 1].sum())))
    support = abs(x - y[..., :x.shape[-1]]).sum(axis=-1) / den
    return dict(scope="Accepted partial or complete matched trajectories, no extrapolation and no relic inference.",
                planned_steps=24, accepted_matched_steps=n, q32=r, q64=s, shared_norms=["number", "energy", "pressure"],
                shared_L1=support.tolist())


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos", type=Path)
    p.add_argument("work", type=Path)
    p.add_argument("steps", type=int)
    a = p.parse_args()
    assert 0 < a.steps <= 24
    print(json.dumps(run(a.eos, a.work, a.steps), indent=2, allow_nan=False))
