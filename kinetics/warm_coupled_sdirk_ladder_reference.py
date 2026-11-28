"""Compare the actual first-, second- and fourth-block endpoints, no fitting."""
import argparse
import json
from pathlib import Path

import numpy as np

from warm_light_occupation_stage_audit_reference import digest,read


def run(root,be):
    locations=[
        (root/"068-abs1",root/"068-abs2",root/"068-abs4"),
        (be/"068-locked1-expand-finite-all",be/"068-locked2-expand-finite-all",root/"068-be4")
    ]
    ref=None
    all_rows=[]
    for names,method in zip(locations,("sdirk2","be"),strict=True):
        rows=[]
        for p,n in zip(names,(1,2,4),strict=True):
            a=read(p/"summary.json")
            v=a["rows"]
            assert len(v)==n
            last=v[-1]
            keys=(a["input_sha256"],a["eos_sha256"],a["cut_GeV"],a["W_initial_over_rho"],
                  last["qsha256"],last["wsha256"],last["c"],last["time_J"])
            if ref is None:
                ref=keys
            assert ref==keys
            for i,r in enumerate(v,1):
                y=read(p/f"step-{i}.json")
                assert y["row"]==r and digest(y["daughter_snapshot"])==r["daughter_hash"]
                assert r.get("collision_clock_method","be")==method
                assert r["solver"]["residual"]<8e-12
            W=np.array([r["W_over_initial_daughter_rho"] for r in v])
            assert np.all(np.diff(np.r_[a["W_initial_over_rho"],W])<0.)
            rows.append(dict(blocks=n,W=float(W[-1]),max_invariant_defect=max(
                    max(abs(r["relative_energy_defect"]),*map(abs,r["relative_N_defects"])) for r in v)))
        w=np.array([r["W"] for r in rows])
        changes=abs(np.diff(w))
        all_rows.append(dict(method=method,grid=rows,absolute_change=changes.tolist(),
                             coarse_to_fine_difference_ratio=float(changes[0]/changes[1])))
    return dict(scope="Observed discrete refinement changes; no asymptotic order or continuum accuracy is asserted.",
                identical_clock_and_chart=ref,methods=all_rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for v in ("root","be"):
        p.add_argument(v,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.root,a.be),indent=2,allow_nan=False))
