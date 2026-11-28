"""Compare actual count effects of frozen operator schedules."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def run(eos, seed, standard, half, sdirk):
    paths=(eos,seed,standard,half,sdirk)
    method=[json.loads(x.read_text()) for x in paths[2:]]
    Y=list(method[0].values())+list(method[1].values())+list(method[2].values())
    labels=[f"BE:{x}" for x in [*method[0],*method[1]]]+[f"SDIRK2:{x}" for x in method[2]]
    Q=np.array([[x["count_relative_to_initial"] for x in m] for m in Y])
    for y in Y:
        assert [v["outer_subintervals"] for v in y]==[1,2,4]
        assert max(max(abs(v[k]) for k in ("energy_ledger","count_ledger")) for v in y)<2e-10
    dif=Q[:,1:]/Q[:,:-1]-1.
    changed=[(Q[j]/Q[i]-1.).tolist() for i,j in ((0,1),(2,1),(3,4),(0,3),(1,4))]
    return dict(methods=labels,n_subintervals=[1,2,4],raw_moment_data=method,
                frozen_pair_and_energy_relative_by_method=Q.tolist(),
                time_grid_relative_by_method=dif.tolist(),
                schedule_change_labels=["BE reverse/standard","BE reverse/two-half-collisions",
                         "SDIRK2 reverse/standard","SDIRK2/BE standard","SDIRK2/BE reverse"],
                scheduling_relative_count_differences=changed,
                sha256={v.name:hashlib.sha256(v.read_bytes()).hexdigest() for v in paths},
                scope="One small normalized interval at fixed mid-grid scale, pressure, visible background and elastic panels. Pair kicks and stiff prompt/elastic maps differ in order and subdivision. Not an expanding-trajectory time test, cost comparison, or proof of general split order.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for x in ("eos","seed","standard","half","sdirk"):
        p.add_argument(x,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.seed,a.standard,a.half,a.sdirk),indent=2,allow_nan=False))
