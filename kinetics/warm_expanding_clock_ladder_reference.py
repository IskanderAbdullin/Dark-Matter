"""Compare actual endpoints of a fixed-parent-chart temporal refinement ladder."""
import argparse
import json
from pathlib import Path

import numpy as np

from warm_light_occupation_stage_audit_reference import digest,read


def analyze(folder,hard,blocks,modes):
    x={}
    rows=[]
    for n in blocks:
        row={}
        for mode in modes:
            p=folder/f"{hard}-locked{n}-{mode}"/"summary.json"
            d=read(p)
            rr=d["rows"]
            assert d["mode"]==mode and len(rr)==n
            assert max(abs(r["relative_energy_defect"]) for r in rr)<3e-10
            for j,r in enumerate(rr,1):
                phase=read(p.parent/f"step-{j}.json")
                assert phase["row"]==r
                assert digest(phase["daughter_snapshot"])==r["daughter_hash"]
                assert r["solver"]["residual"]<8e-12
            a=rr[-1]
            invariants=(a["qsha256"],a["wsha256"],a["c"],d["input_sha256"],a["time_J"],
                        d["W_initial_over_rho"],d["cut_GeV"])
            if not x:
                x=dict(data=invariants,W_initial=d["W_initial_over_rho"])
            assert invariants==x["data"]
            row[mode]=dict(W=a["W_over_initial_daughter_rho"],e=a["relative_energy_defect"],
                           peak_step_residual=max(r["solver"]["residual"] for r in rr),
                           daughter_hash=a["daughter_hash"],
                           decreasing_at_every_sample=bool(np.all(np.diff(
                               [d["W_initial_over_rho"]]+[r["W_over_initial_daughter_rho"] for r in rr])<0.)))
        rows.append(dict(outer_blocks=n,ends=row))
    for mode in modes:
        w=np.array([r["ends"][mode]["W"] for r in rows])
        assert np.all((w>0.)&(w<x["W_initial"]))
        for r,old,new in zip(rows[1:],w[:-1],w[1:],strict=True):
            r["ends"][mode]["relative_shift_from_preceding_grid"]=float(new/old-1.)
    if len(modes)==2:
        for r in rows:
            a,b=[r["ends"][m]["W"] for m in modes]
            r["second_vs_first_mode_relative_shift"]=b/a-1.
    return dict(scope="Matched numerical endpoints and sampled changes only; no extrapolated continuum value or accuracy enclosure.",
                clock_and_chart=x,rows=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("folder",type=Path)
    p.add_argument("hard")
    p.add_argument("blocks",type=int,nargs="+")
    p.add_argument("--finite-only",action="store_true")
    a=p.parse_args()
    modes=("expand-finite-all",) if a.finite_only else ("expand-stationary-all","expand-finite-all")
    print(json.dumps(analyze(a.folder,a.hard,a.blocks,modes),indent=2,allow_nan=False))
