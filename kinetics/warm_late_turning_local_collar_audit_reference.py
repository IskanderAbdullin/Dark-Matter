"""Check hashes, matched initial moments and ledger for the local collar replay."""

import argparse
import json
from pathlib import Path

import numpy as np

from late_radiation_background_reference import TabulatedRadiation
from warm_coupled_serial_restart_reference import sha
from warm_late_turning_local_collar_replay_reference import remap
from warm_light_occupation_stage_audit_reference import digest, read
from warm_prompt_reentry_shape_audit_reference import get_mom


def audit(eos, original, refined):
    raw = [read(original / f"step-{i}.json") for i in range(7)]
    fine = [read(refined / f"step-{i}.json") for i in range(1, 7)]
    short = read(refined / "report.json")
    for step in raw + fine:
        assert digest(step["snapshot"]) == step["hash"]
        assert step["plan"]["eos_sha256"] == sha(eos)
    assert short["rows"] == [{k:v for k,v in step.items() if k != "snapshot"} for step in fine]
    initial = raw[0]["snapshot"]
    p = np.array(initial["lab_momentum_edges_GeV"])
    q = np.array(fine[0]["snapshot"]["lab_momentum_edges_GeV"])
    assert all(np.any(q == x) for x in p)
    split = q[~np.isin(q, p)]
    assert len(split) == 4 and np.all((split >= .645) & (split <= .765))
    def exact(step, edges):
        s = step["snapshot"]
        assert np.array_equal(edges, s["lab_momentum_edges_GeV"])
        return get_mom(edges, [np.array(s["entropy_variable_"+k]) for k in ("scalar", "fermion")])
    n, _, e = exact(raw[0], p)
    N, _, E = get_mom(q, remap(p, q, initial))
    # Each refined shell is contained within one and only one original shell.
    i = np.searchsorted(p, (q[:-1]+q[1:])/2, side="right")-1
    max_restrict = 0.
    for a,b in ((n,N), (e,E)):
        for m0,m1 in zip(a,b):
            reconcile = np.bincount(i, weights=m1, minlength=len(m0))
            assert np.allclose(reconcile,m0,rtol=2e-12,atol=0.)
            max_restrict = max(max_restrict, float(np.sum(abs(reconcile-m0))/sum(m0)))
    rad = TabulatedRadiation(eos,initial["T_visible_GeV"]).state(1.)[1]
    dh = abs(np.sqrt((rad+E.sum())/(rad+e.sum()))-1.)
    assert dh < 2e-14
    rows=[]
    for a,b in zip(raw[1:],fine):
        assert b["guide_was_initial_guess_sha256"] == a["hash"]
        assert b["plan"]["input_sha256"] == raw[0]["hash"]
        assert b["plan"]["cut_GeV"] == raw[0]["plan"]["normalized_cut_GeV"]
        assert b["time"] == a["time"] and b["same_clock_original_tail_excess"] == a["tail_energy"]
        n0,_,e0=exact(a,p)
        n1,_,e1=exact(b,q)
        error=[]
        for nn,ee in ((n0,e0),(n1,e1)):
            error.append([float(max(abs(nn.sum(axis=1)/n.sum(axis=1)-1.))),
                          float(abs(ee.sum()/e.sum()-1.))])
        assert np.max(error) < 2e-12
        assert np.sign(a["signed_native_flux"]) == np.sign(b["native_signed_flux"])
        rows.append(dict(time=b["time"], original_and_split_W=[a["tail_energy"],b["tail_excess_over_rho"]],
                         original_and_split_flux=[a["signed_native_flux"],b["native_signed_flux"]],
                         original_and_split_number_and_energy_defects=error))
    bracket=[(a["time"],b["time"]) for a,b in zip(rows,rows[1:])
             if a["original_and_split_flux"][0]*b["original_and_split_flux"][0]<0.]
    assert len(bracket)==1
    peak=max(row["original_and_split_W"][0] for row in rows)
    return dict(original_initial_hash=raw[0]["hash"], split_points=split.tolist(),
                relative_H_difference=float(dh), largest_shell_restriction_L1=max_restrict,
                common_flux_sign_bracket=bracket, max_W_shift_over_original_peak=max(
                abs(a-b)/peak for a,b in (row["original_and_split_W"] for row in rows)), steps=rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("eos","original","refined"):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    print(json.dumps(audit(a.eos,a.original,a.refined),indent=2,allow_nan=False))
