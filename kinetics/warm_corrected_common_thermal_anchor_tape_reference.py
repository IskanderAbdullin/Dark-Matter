"""Checkpoint common-equilibrium tangents at fixed corrected trajectory stages."""

import argparse
import json
from pathlib import Path

from warm_coupled_serial_restart_reference import atomic, sha
from warm_light_occupation_stage_audit_reference import digest, read
from warm_pair_common_moment_equilibrium_reference import report


def run(eos,wide,out):
    paths={0:wide/"guide-start.json",**{i:wide/"q64"/f"step-{i}.json" for i in (6,12,18,24)}}
    fingerprint=dict(eos_sha256=sha(eos),module_sha256=sha(Path(report.__code__.co_filename)),
                     relaxation_sha256=sha(Path(__file__).with_name("warm_common_pair_thermal_relaxation_reference.py")))
    out.mkdir(parents=True,exist_ok=True)
    rows=[]
    for i,p in paths.items():
        src=read(p)
        h=src["snapshot_sha256"]
        assert digest(src["snapshot"])==h
        target=out/f"fit-{i}.json"
        if not target.exists():
            v=report(eos,p,True)
            assert v["snapshot_sha256"]==h
            atomic(target,dict(plan=fingerprint,result=v))
        row=read(target)
        assert row["plan"]==fingerprint and row["result"]["snapshot_sha256"]==h
        rows.append(row["result"])
    return dict(scope="Independently fitted common moments; no trajectory is assumed thermal.",
                plan=fingerprint,cases=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for x in ("eos","wide","out"):
        p.add_argument(x,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.wide,a.out),indent=2,allow_nan=False))
