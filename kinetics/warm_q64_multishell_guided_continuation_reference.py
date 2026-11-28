"""Guide the wider nine-cell trajectory using accepted corrected q32 seeds only."""

import argparse
import json
from pathlib import Path

import coupled_warm_pair_continuation_reference as kernel
import warm_cold_light_pair_partition_reference as pair
from warm_coupled_serial_restart_reference import atomic, sha
from warm_light_occupation_stage_audit_reference import digest, read
from warm_q64_guided_serial_reference import run


def execute(eos, warm64, narrow, path, count):
    path.mkdir(parents=True, exist_ok=True)
    seed = read(warm64)[-1]
    assert seed["snapshot_sha256"] == digest(seed["snapshot"])
    guide = [read(narrow / f"step-{i}.json") for i in range(1,count+1)]
    start = path / "guide-start.json"
    block = path / "guide-profile.json"
    for f, data in ((start,seed),(block,dict(blocks=[{"checkpoint":g} for g in guide]))):
        if f.exists():
            assert read(f) == data
        else:
            atomic(f,data)
    method = Path(pair.__file__)
    tag_order = (4,8,6,8,8)
    plan = path / "correction-manifest.json"
    row = dict(cut_cells=9,tagged_order=list(tag_order), input_sha256=sha(warm64), eos_sha256=sha(eos),
               guide_sha256=sha(block), corrections_sha256={j.name:sha(j) for j in
                  [method,method.with_name("warm_cold_light_finite_pair_reference.py"),
                   method.with_name("archived_scalar_pair_reaction_reference.py"),Path(__file__)]})
    if plan.exists():
        assert read(plan)==row
    else:
        assert not list((path/"q64").glob("step-*.json"))
        atomic(plan,row)
    old = kernel.rate

    def rate(p, th, o):
        V, e, ratio, _ = pair.split_rate(p,th,o,tag_order,cut_cells=9)
        return V,e,ratio

    try:
        kernel.rate=rate
        data=run(eos,start,block,path/"q64")
    finally:
        kernel.rate=old
    return dict(scope="Independent q64 targets; corrected q32 states are nonlinear initial guesses only.",
                correction_manifest=row, continuation=data)


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for j in ("eos","warm64","narrow","path"):
        p.add_argument(j,type=Path)
    p.add_argument("--blocks",type=int,required=True)
    a=p.parse_args()
    print(json.dumps(execute(a.eos,a.warm64,a.narrow,a.path,a.blocks),indent=2,allow_nan=False))
