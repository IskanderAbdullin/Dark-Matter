"""Compare actual early-refined and uniform accepted clocks at one fixed cost."""

import argparse
import json
from pathlib import Path

import numpy as np

from warm_clock_path_error_profile_reference import audited


def run(uniform12,early12,uniform24):
    groups=[uniform12,early12,uniform24]
    solved=[audited(p) for p in groups]
    keys=[(r["input_sha256"],r["eos_sha256"],r["W_initial_over_rho"],r["mode"],r["cut_GeV"],
           *(r["rows"][-1][k] for k in ("qsha256","wsha256","c","time_J"))) for r in solved]
    assert len(set(keys))==1
    assert [len(r["rows"]) for r in solved]==[12,12,24]
    anchor=solved[0]["W_initial_over_rho"]
    ends=[]
    for folder,s in zip(groups,solved,strict=True):
        W=np.array([x["W_over_initial_daughter_rho"] for x in s["rows"]])
        h=np.array([x["step_H"] for x in s["rows"]])
        assert np.all(h>0.) and W[-1]>0.
        ends.append(dict(case=folder.name,W=float(W[-1]),normalized_clock_steps=(h/h.sum()).tolist(),
                         decline_every_sample=bool(np.all(np.diff(np.r_[anchor,W])<0.))))
    same=[]
    for j in range(1,5):
        p=f"step-{j}.json"
        same.append((early12/p).read_bytes()==(uniform24/p).read_bytes())
    assert all(same)
    x,y,z=[p["W"] for p in ends]
    return dict(scope="Measured difference to a finer positive clock, not a certified continuum error.",
                initial_W=anchor,inputs=keys[0],cases=ends,shared_first_four_blocks_byte_identical=same,
                early_vs_uniform12_relative=y/x-1.,
                relative_gaps_to_uniform24=[x/z-1.,y/z-1.],
                early_to_uniform12_absolute_gap_ratio=abs(y-z)/abs(x-z))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("uniform12","early12","uniform24"):
        p.add_argument(name,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.uniform12,a.early12,a.uniform24),indent=2,allow_nan=False))
