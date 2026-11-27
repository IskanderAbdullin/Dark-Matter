"""Separate endpoint timestep error from the error along matched saved phases."""

import argparse
import json
from pathlib import Path

import numpy as np

from warm_light_occupation_stage_audit_reference import digest,read


def audited(root):
    s=read(root/"summary.json")
    for j,r in enumerate(s["rows"],1):
        raw=read(root/f"step-{j}.json")
        assert raw["row"]==r and r["daughter_hash"]==digest(raw["daughter_snapshot"])
        assert abs(r["relative_energy_defect"])<3e-10
    return s


def run(a,b):
    s,t=audited(a),audited(b)
    for k in ("input_sha256","eos_sha256","W_initial_over_rho","mode","cut_GeV"):
        assert s[k]==t[k]
    r0,r1=s["rows"],t["rows"]
    for k in ("qsha256","wsha256","c","time_J"):
        assert r0[-1][k]==r1[-1][k]
    rows=[]
    w=s["W_initial_over_rho"]
    for x in r0:
        i=int(np.argmin([abs(r["time_J"]-x["time_J"]) for r in r1]))
        y=r1[i]
        if abs(y["time_J"]-x["time_J"])>3e-21:
            continue
        v,u=x["W_over_initial_daughter_rho"],y["W_over_initial_daughter_rho"]
        assert min(u,v)>0.
        rows.append(dict(time_H=x["time_J"],coarse_W=v,refined_W=u,
                         difference_over_initial_W=(v-u)/w,
                         difference_over_matched_refined_W=v/u-1.))
    assert rows and rows[-1]["time_H"]==r0[-1]["time_J"]
    return dict(scope="Signed differences on recorded endpoints only, not a local error estimator or an error enclosure.",
                starting_daughter_hash=s["input_sha256"],W_initial_over_rho=w,
                coarse_and_refined_counts=[len(r0),len(r1)],samples=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("coarse",type=Path)
    p.add_argument("refined",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.coarse,a.refined),indent=2,allow_nan=False))
