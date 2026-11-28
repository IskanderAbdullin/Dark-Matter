"""Separate elastic, lifted pair and prompt rates on accepted cold-corrected states."""

import argparse
import json
from pathlib import Path

import warm_recoil_endpoint_channel_budget_reference as budget
from warm_cold_light_pair_partition_reference import split_rate
from warm_light_occupation_stage_audit_reference import digest, read


def measure(eos, path):
    row = read(path)
    if isinstance(row,list):
        row = row[-1]
    S = row["snapshot"]
    assert row["snapshot_sha256"] == digest(S)
    base = budget.rate

    def corrected(p, th, order):
        V, err, ratio, _ = split_rate(p,th,order,(4,8,6,8,8),cut_cells=9)
        return V,err,ratio

    try:
        budget.rate = corrected
        a = budget.at(eos,S,[((3,2,2,8),(3,7,6,2)),((6,4,4,12),(4,10,10,4))])
    finally:
        budget.rate = base
    assert all(x["sum_per_log_a"] < 0 for x in a)
    return dict(T_GeV=S["T_visible_GeV"],snapshot_sha256=digest(S),from_archive=path.name,cases=a)


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("snapshots",nargs="+",type=Path)
    a=p.parse_args()
    print(json.dumps([measure(a.eos,f) for f in a.snapshots],indent=2,allow_nan=False))
