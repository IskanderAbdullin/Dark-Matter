"""Replay the common native cold spike at fixed inputs using the split gain."""
import argparse
import json
from pathlib import Path
import numpy as np
import coupled_warm_pair_continuation_reference as kernel
from archived_scalar_pair_reaction_reference import __file__ as gross_path
from warm_cold_light_pair_partition_reference import split_rate
from warm_light_occupation_stage_audit_reference import read, digest
from warm_coupled_serial_restart_reference import atomic, sha
from warm_history_lower_matched_audit_reference import measures
from warm_q64_history_serial_reference import run

p=argparse.ArgumentParser(description=__doc__)
p.add_argument("eos", type=Path)
p.add_argument("work", type=Path)
p.add_argument("out", type=Path)
p.add_argument("--q", choices=(32,64), type=int, required=True)
p.add_argument("--fine", action="store_true")
p.add_argument("--subdivisions", choices=(1,2,4), type=int, default=1)
a=p.parse_args()
a.out.mkdir(parents=True, exist_ok=True)
A,B,C=[read(a.work/f"q{a.q}"/f"step-{i}.json") for i in (10,11,12)]
S=B["snapshot"]
assert A["snapshot_sha256"]==B["previous"] and B["snapshot_sha256"]==C["previous"]
seed=a.out/"seed.json"
atomic(seed,[A,B])
b=C["plan"]["step_expansion"]**(1/a.subdivisions)
assert sha(a.eos) == C["plan"]["EOS_sha256"]
mode=(6,12,9,12,12) if a.fine else (4,8,6,8,8)
logs=[]
original=kernel.rate

def callback(p,th,o):
    rr,ledger,ratio,meta=split_rate(p,th,o,mode)
    logs.append(meta)
    return rr,ledger,ratio

try:
    report=run(a.eos,seed,b,a.subdivisions,a.out/"native")
    N=read(a.out/"native"/f"step-{a.subdivisions}.json")
    assert report["final_snapshot_sha256"] == N["snapshot_sha256"]
    if a.subdivisions == 1:
        assert N["snapshot_sha256"] == C["snapshot_sha256"]
    kernel.rate=callback
    corrected=run(a.eos,seed,b,a.subdivisions,a.out/"lifted")
finally:
    kernel.rate=original
R=read(a.out/"lifted"/f"step-{a.subdivisions}.json")
assert corrected["final_snapshot_sha256"] == R["snapshot_sha256"]
assert abs(R["snapshot"]["T_visible_GeV"]/C["snapshot"]["T_visible_GeV"] - 1.) < 3e-14
def occ(x):
    e=np.asarray(x["lab_momentum_edges_GeV"])
    return measures(x)[0,1,0]/(4*e[1]**3/(6*np.pi**2))
def row(x):
    s=x["snapshot"]
    assert digest(s)==x["snapshot_sha256"]
    assert abs(x["operator_summary"]["count_ledger"]) < 2e-11
    return dict(T_GeV=s["T_visible_GeV"],sha256=x["snapshot_sha256"],occupation=float(occ(s)),
                slope=float((occ(s)-occ(S))/(a.subdivisions*np.log(b))),
                energy_defect=x["operator_summary"]["energy_ledger"],
                total_count_defect=x["operator_summary"]["count_ledger"])
from warm_cold_light_pair_partition_reference import __file__ as method
delta = measures(R["snapshot"]) - measures(N["snapshot"])
scale = measures(N["snapshot"])
shapes = (np.abs(delta).sum(axis=-1) / scale.sum(axis=-1)).tolist()
n = row(N)
r = row(R)
if a.subdivisions == 1:
    assert n["slope"] > 0
assert r["slope"] < 0
print(json.dumps(dict(q=a.q, subdivisions=a.subdivisions, starting_sha256=B["snapshot_sha256"], lifts=logs,
    order=mode, lift_sha256=sha(Path(method)), gross_sha256=sha(Path(gross_path)),
    normalized_shell_L1_number_energy_pressure=shapes, native=n, corrected=r),indent=2,allow_nan=False))
