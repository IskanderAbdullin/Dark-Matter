"""Verify matched locked stiff clocks and distinguish fallback from pair order."""

import argparse
import json
from pathlib import Path

import numpy as np

from warm_history_lower_matched_audit_reference import measures
from warm_light_occupation_stage_audit_reference import digest,read


def run(nine,free_default,free_fine,locked_default,locked_fine):
    R=[read(p) for p in (nine,free_default,free_fine,locked_default,locked_fine)]
    S=[r["snapshot"] for r in R]
    for r,s in zip(R,S,strict=True): assert r["snapshot_sha256"]==digest(s)
    A,B,C,D,E=R
    b,d=D["source"],E["source"]
    assert b["parent_sha256"]==d["parent_sha256"]==A["parent_snapshot_sha256"]
    assert b["nine_cell_sha256"]==d["nine_cell_sha256"]==A["snapshot_sha256"]
    assert b["guide_sha256"]==B["snapshot_sha256"] and d["guide_sha256"]==C["snapshot_sha256"]
    assert b["chart_sha256"]==d["chart_sha256"] and b["eos_sha256"]==d["eos_sha256"]
    assert b["order"]==[3,7,6,2] and d["order"]==[4,10,10,4]
    steps=[v["logged_stiff_substeps"] for v in (D,E)]
    assert len(steps[0])==len(steps[1])==b["selected_substeps"]==d["selected_substeps"]
    clocks=[[s["step_H"] for s in j] for j in steps]
    assert clocks[0]==clocks[1] and all(h==clocks[0][0] for h in clocks[0])
    for r,s in zip(R[1:],S[1:],strict=True):
        np.testing.assert_allclose(s["T_visible_GeV"],S[0]["T_visible_GeV"],rtol=3e-14)
        np.testing.assert_allclose(s["lab_momentum_edges_GeV"],S[0]["lab_momentum_edges_GeV"],rtol=3e-14)
    X=[measures(r) for r in S]
    n=X[0].sum(axis=-1)
    L=lambda x,y:(abs(X[x]-X[y]).sum(axis=-1)/n).tolist()
    return dict(scope="Fixed equal implicit microsteps for both charts and identical outer clock. "
                      "Selected step only, not a whole-trajectory convergence proof.",
                snapshot_sha256=[r["snapshot_sha256"] for r in R],
                metric_order=["number","full_energy","pressure"],substep_sizes_in_Hubble_units=clocks[0],
                original_adaptive_spread=L(1,2),matched_locked_spread=L(3,4),
                default_locked_vs_adaptive=L(3,1),fine_locked_vs_adaptive=L(4,2),
                locked_cold_slopes=[x["cold_slope_per_ln_a"] for x in (D,E)],
                residuals=[[z["residual"] for z in s] for s in steps],
                budgets=[r["energy_and_count_ledger"] for r in (D,E)])


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("nine","free_default","free_fine","locked_default","locked_fine"):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.nine,a.free_default,a.free_fine,a.locked_default,a.locked_fine),indent=2,allow_nan=False))
