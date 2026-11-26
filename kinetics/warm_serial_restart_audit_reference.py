"""Compare checkpointed warm output to an unbroken matched outer grid."""
import argparse
import json
from pathlib import Path

import numpy as np

from warm_coupled_serial_restart_reference import bytes_of, sha
from warm_exterior_restart_audit_reference import ext, public_count, shell


def read(p):
    return json.loads(p.read_text())


def compare(start, first, finished, broken, unbroken):
    files=[start,first,finished,broken,unbroken]
    S,A,B,C,D=map(read,files)
    assert A["plan"]==B["plan"] and A["plan_sha256"]==B["plan_sha256"]
    assert A["blocks_finished"]==1 and B["blocks_finished"]==B["blocks_total"]==2
    assert A["ledger"]==B["ledger"][:1] and np.max([abs(q["energy_ledger"]) for q in B["ledger"]])<2e-10
    assert bytes_of(B["snapshot"])==bytes_of(C["snapshots"][-1])
    a,b=map(lambda x:x["snapshots"][-1],(C,D))
    x,y=map(shell,(a,b))
    err=float(np.max(np.sum(abs(x-y),axis=2)/np.sum(abs(y),axis=2)))
    E,F=map(ext,(a,b))
    en=float(abs(E.energy(1.)/F.energy(1.)-1.))
    tr=float(abs(E.trace/F.trace-1.))
    u=float(abs(E.upper_bound(1.)/F.upper_bound(1.)-1.))
    z=float(abs(a["T_visible_GeV"]/b["T_visible_GeV"]-1.))
    initial=public_count(S["snapshots"][-1])
    count=float(max(abs(public_count(v)/initial-1.) for v in (a,b)))
    assert max(err,en,tr,u,z,count)<2e-10
    return dict(resumed_blocks=B["blocks_finished"],final_T_GeV=a["T_visible_GeV"],
                energy_ledger_per_resumed_block=[q["energy_ledger"] for q in B["ledger"]],
                relative_shell_count_energy_defect=err,relative_exterior_mean_energy_defect=en,
                relative_exterior_upper_energy_defect=u,relative_exterior_trace_defect=tr,
                relative_clock_defect=z,maximum_global_count_defect=count,
                same_plan_sha256=B["plan_sha256"],record_sha256={p.name:sha(p) for p in files},
                scope="Two-block restart experiment on a short prescribed warm grid, not a new long warm trajectory or convergence test.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for x in ("start","first","finished","serial","unbroken"):
        p.add_argument(x,type=Path)
    a=p.parse_args()
    print(json.dumps(compare(a.start,a.first,a.finished,a.serial,a.unbroken),indent=2,allow_nan=False))
