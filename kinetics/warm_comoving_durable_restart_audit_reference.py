"""Check byte-identical resumption and refusal of inconsistent saved records."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import tempfile

import numpy as np

from warm_comoving_durable_restart_reference import restore
from warm_coupled_serial_restart_reference import atomic
from warm_light_occupation_stage_audit_reference import read


def h(v):
    return hashlib.sha256(np.asarray(v,dtype="<f8").tobytes()).hexdigest()


def audit(reference,resumed):
    p=resumed/"restart-state.json"
    source=read(p)["payload"]
    n=source["accepted_blocks"]
    r=read(resumed/"summary.json")
    assert n==len(r["rows"])
    names=[f"step-{j}.json" for j in range(1,n+1)]+["summary.json"]
    identical=[]
    for name in names:
        a=(reference/name).read_bytes()
        b=(resumed/name).read_bytes()
        assert a==b
        identical.append(dict(file=name,sha256=hashlib.sha256(a).hexdigest(),bytes=len(a)))
    scales=np.r_[1.,[x["c"] for x in r["rows"]]]
    lineage=source["key"]
    hv=source["parent_occupations"]
    Q=np.empty_like(np.asarray(hv))
    assert restore(p,lineage,scales,r["input_sha256"],Q,h)[0]==source
    try:
        restore(p,"other-lineage",scales,r["input_sha256"],Q,h)
    except AssertionError:
        different_rejected=True
    else:
        different_rejected=False
    assert different_rejected
    with tempfile.TemporaryDirectory() as t:
        box=Path(t)
        for name in names+["restart-state.json"]:
            shutil.copyfile(resumed/name,box/name)
        broken=read(box/"step-1.json")
        broken["row"]["W_over_initial_daughter_rho"]*=1.000001
        atomic(box/"step-1.json",broken)
        before=(box/"step-1.json").read_bytes()
        try:
            restore(box/p.name,lineage,scales,r["input_sha256"],Q,h)
        except AssertionError:
            corrupted_rejected=True
        else:
            corrupted_rejected=False
        assert corrupted_rejected and before==(box/"step-1.json").read_bytes()
    return dict(scope="Byte identity tested on a single locked cloud CPU run, not portability across library versions.",
                initial_hash=r["input_sha256"],accepted_blocks=n,parent_nodes=len(hv),
                checkpoint_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
                accepted_exterior_cohorts=len(source["exterior"]["cohorts"]),outputs=identical,
                changed_lineage_rejected=different_rejected,changed_previous_output_rejected=corrupted_rejected)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("reference",type=Path)
    p.add_argument("resumed",type=Path)
    a=p.parse_args()
    print(json.dumps(audit(a.reference,a.resumed),indent=2,allow_nan=False))
