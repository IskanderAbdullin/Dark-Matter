"""Conservative hybrid scalar-pair projection with tagged inverse cold events."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from archived_scalar_pair_reaction_reference import gross
from joint_prompt_two_moment_evolution_reference import Cells
from joint_visible_heavy_adiabatic_reference import snapshot_functions
from scalar_pair_expanding_continuation_reference import e_total
from three_species_decay_onsager_reference import m, mp
from warm_cold_light_finite_pair_reference import tagged
from warm_light_occupation_stage_audit_reference import digest, read


def lifted(p, F, B, fe, be, order):
    """Each event has at least one outgoing fermion in the first cell.

    Tag either charge and spin, with four possible tags in the radial
    measure. Exactly two tags correspond to the same event if both light
    legs lie in the cold cell; the division below restores event counting.
    """
    a, r, i, o, z = order
    ecut = fe[1]
    cold = np.zeros((2, 2, len(p) + 1))
    mass = (mp, m)

    def scat(j, E, w):
        E = np.broadcast_to(E, w.shape).ravel()
        w = w.ravel()
        if not np.any(w):
            return
        ix = np.searchsorted((be, fe)[j], E, side="right")
        cold[j, 0] += np.bincount(ix, w, minlength=len(p) + 1)
        cold[j, 1] += np.bincount(ix, w * np.maximum(0., E - mass[j]), minlength=len(p) + 1)

    x, w = roots_legendre(a)
    q = p[1] * (1 + x) / 2
    weights = 4 * p[1] / 2 * q * q * w / (2 * np.pi**2)
    expected = 0.
    both = 0.
    for pt, volume in zip(q, weights, strict=True):
        def emit(E, F, e0, e1, Q):
            nonlocal both
            event = volume * Q / (1 + (F < ecut))
            both += np.sum(event * (F < ecut))
            scat(0, e0, event)
            scat(0, e1, event)
            scat(1, E, event)
            scat(1, F, event)

        g, _ = tagged(pt, F, B, fe, be, (r, i, o, z), emit=emit)
        expected += volume * g[0]
    check = abs(cold[1, 0, 1] / expected - 1.)
    assert check < 3e-13
    # Separately close the projected channel before adding its complement.
    mag = cold[0, 0].sum()
    diff = cold * np.array([-1, 1])[:, None, None]
    assert abs(diff[:, 0].sum() / mag) < 4e-13
    assert abs(e_total(diff) / (mag * m)) < 4e-13
    return cold, float(both), float(check)


def split_rate(p, theta, native_order, tag_order=(4, 8, 6, 8, 8), reverse=True):
    B, be, F, fe = snapshot_functions(dict(lab_momentum_edges_GeV=p,
        entropy_variable_scalar=theta[0], entropy_variable_fermion=theta[1]))
    _, Fg, L = gross(p, m, mp, F, B, False, *native_order, project=True)
    if reverse:
        _, Rg, C = gross(p, mp, m, B, F, True, *native_order, project=True,
                        min_final_energy=fe[1])
        cold, both, check = lifted(p, F, B, fe, be, tag_order)
        assert not np.any(C[1, :, 1])
    else:
        C, cold = np.zeros_like(L), np.zeros_like(L)
        Rg = np.zeros_like(Fg)
        both = check = 0.
    G = C + cold
    net = np.array([L[1] - G[0], G[1] - L[0]])
    nscale = 2 * (Fg[1] + Rg[1]) + cold[0, 0].sum()
    escale = Fg[2] + Rg[2] + cold[0,1].sum() + mp*cold[0,0].sum()
    dN = abs(net[:, 0].sum()) / max(nscale, 1e-280)
    dE = abs(e_total(net)) / max(escale, 1e-280)
    assert max(dN, dE) < 2e-11 and np.all(net[:, :, [0,-1]] >= 0.)
    events = Rg[1] + cold[0, 0].sum()/2
    info = dict(total_number_ledger=float(dN), total_energy_ledger=float(dE), tag_closure=check,
                cold_reverse_event_fraction=float(cold[0,0].sum()/2/max(events,1e-280)),
                both_cold_event_fraction=float(both/max(events,1e-280)),
                cold_gain=float(cold[1,0,1]), cold_loss=float(L[0,0,1]))
    return net, max(dN,dE), events/Fg[1], info


def control(files):
    rows = []
    for file in files:
        row = read(file)
        if isinstance(row, list):
            row = row[-1]
        S = row["snapshot"]
        assert row["snapshot_sha256"] == digest(S)
        p = np.asarray(S["lab_momentum_edges_GeV"])
        th = [S["entropy_variable_scalar"], S["entropy_variable_fermion"]]
        v = Cells(p,m,4).vol[0]
        choices = []
        for o in ((4,8,6,8,8),(6,12,9,12,12)):
            n,_,_,d = split_rate(p, th, (3,7,6,2), o)
            choices.append(dict(tag_order=o, net_cold_rate_GeV=float(n[1,0,1]/v), **d))
        assert abs(choices[-1]["cold_gain"]/choices[0]["cold_gain"] - 1.) < 1.e-6
        rows.append(dict(T_GeV=S["T_visible_GeV"], snapshot_sha256=digest(S), projection=choices))
    return dict(scope="Frozen hybrid collision and four-leg moment conservation; not evolved.", cases=rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("files", nargs="+", type=Path)
    print(json.dumps(control(p.parse_args().files), indent=2, allow_nan=False))
