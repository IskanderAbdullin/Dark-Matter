"""Compare two completed moment-matched hard cohort realizations."""

import argparse
import json
from pathlib import Path

import numpy as np

from frw_exterior_restart_reference import import_restart
from three_species_decay_onsager_reference import m, mp
from warm_history_lower_matched_audit_reference import measures
from warm_light_occupation_stage_audit_reference import digest, read


def compare(paths, metrics):
    rows = [read(p) for p in paths]
    a, b, vac = S = [r["snapshot"] for r in rows]
    for r,s in zip(rows,S,strict=True):
        assert r["snapshot_sha256"]==digest(s)
        assert r["plan"]["outer_log_scale_grid"]==rows[0]["plan"]["outer_log_scale_grid"]
        assert r["plan"]["eos_sha256"]==rows[0]["plan"]["eos_sha256"]
        np.testing.assert_allclose(s["T_visible_GeV"],vac["T_visible_GeV"],rtol=2e-14)
        np.testing.assert_allclose(s["scale_factor"],vac["scale_factor"],rtol=2e-14)
        np.testing.assert_allclose(s["lab_momentum_edges_GeV"],vac["lab_momentum_edges_GeV"],rtol=2e-14)
    V = measures(vac)
    A, B = measures(a), measures(b)
    D = abs(A-B).sum(axis=-1)/V.sum(axis=-1)
    a1,b1 = [import_restart(s["positive_exterior_restart"]).array() for s in (a,b)]
    added = np.array([[q[5,q[1]==d].sum() for d in (mp,m)] for q in (a1,b1)])/V[0].sum(axis=-1)
    M = [read(x) for x in metrics]
    for u,r,hard in zip(M,rows[:2],(.52,.58),strict=True):
        assert u["hard_GeV"]==hard
        assert u["results"]["terminal_q64_sha256"]==r["snapshot_sha256"]
        assert u["results"]["vacuum_q64_sha256"]==rows[-1]["snapshot_sha256"]
    def rate(k):
        return np.array([u["results"][k][-1] for u in M])
    k="high_inverse_and_tagged_energy_rate"
    I, F = rate("initial_"+k), rate("terminal_"+k)
    return dict(scope="Two prescribed realizations of the same saved cohorts, not a tail bound.",
                metric_order=["number","full_energy","pressure"],
                first_stage_hardness=[.52,.58],
                relative_L1_same_grid=D.tolist(), new_receipt_count_relative_to_vacuum=added.tolist(),
                initial_ratio_hard052_over_hard058=(I[0]/I[1]).tolist(),
                terminal_ratio_hard052_over_hard058=(F[0]/F[1]).tolist(),
                terminal_over_initial=(F/I).tolist(),
                terminal_sha256=[r["snapshot_sha256"] for r in rows])


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("terminal052","terminal058","vacuum","metrics052","metrics058"):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    print(json.dumps(compare((a.terminal052,a.terminal058,a.vacuum),(a.metrics052,a.metrics058)),indent=2,allow_nan=False))
