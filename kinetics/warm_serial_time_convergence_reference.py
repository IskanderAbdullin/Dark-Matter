"""Audit equal-support time refinement on the archived warm short benchmark."""
import argparse
import json
import gzip
from pathlib import Path

import numpy as np

from warm_early_serial_support_audit_reference import read, rows, sha


def run(archive, seed, coarse_support, checkpoints):
    grid_sizes=[4*2**j for j in range(len(checkpoints))]
    s=read(seed)["snapshots"][-1]
    a=read(archive)["cases"]
    selected=[v for v in a if v["final_visible_temperature_GeV"]==s["T_visible_GeV"]]
    assert len(selected)==1 and len(checkpoints)>=3
    C,=selected
    src=sha(seed)
    Q=[json.loads(gzip.decompress(p.read_bytes())) if p.suffix=='.gz' else read(p) for p in checkpoints]
    p=Q[0]["plan"]
    for j,(n,z) in enumerate(zip(grid_sizes,Q,strict=True)):
        assert z["blocks_total"]==z["blocks_finished"]==n
        A=z["plan"]
        for key in ("source_sha256","eos_sha256","local_prompt_jacobian","reference_sha256","schema"):
            assert A[key]==p[key]
        assert np.max(abs(np.asarray(A["outer_log_scale_grid"])[::2**j]-p["outer_log_scale_grid"]))<2e-15
    V=[rows(s,z,C,src)[2**j-1::2**j] for j,z in enumerate(Q)]
    assert all(len(x)==4 for x in V)
    ref=read(coarse_support)["matched"]
    result=[]
    for z,j in zip(zip(*V),ref,strict=True):
        b=z[0]["b"]
        T=z[0]["T_GeV"]
        assert b==j["b"] and abs(T/j["T_GeV"]-1.)<1e-12
        assert max(abs(x["b"]/b-1.) for x in z)<1e-13 and max(abs(x["T_GeV"]/T-1.) for x in z)<1e-12
        B=np.asarray([x["inclusive_fixed_count_omega"] for x in z])
        assert np.max(abs(B[0]/j["q16"]["inclusive_fixed_count_omega"]-1.))<1e-12
        D=np.diff(B,axis=0)
        A=np.log2(abs(D[:-1]/D[1:]))
        rel=B[1:]/B[:-1]-1.
        support=np.asarray(j["relative_q32_minus_q16_inclusive"])
        result.append(dict(b=b,T_GeV=T,inclusive_omegas_by_grid=B.tolist(),physical_delta_over_previous_grid=rel.tolist(),
                          observed_order_triples=A.tolist(),same_signed_refinement=(D[:-1]*D[1:]>0.).tolist(),
                          coarse_support_q32_vs_q16=support.tolist(),finest_time_change_over_coarse_support=abs(rel[-1]/support).tolist()))
    return dict(block_counts=grid_sizes,matched=result,
         checkpoint_sha256={p.name:sha(p) for p in checkpoints},recipe_sha256=[sha(p) for p in (archive,seed,coarse_support)],
         scope="Same q16 momentum support, signed operators, initial state and clock; matched prescribed steps and terminal counts only. Empirical apparent order, not a continuum bound or complete phase-space/time convergence, and not a q32 finer-grid comparison.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for x in ("archive","seed","coarse_support"):
        p.add_argument(x,type=Path)
    p.add_argument("checkpoints",type=Path,nargs="+")
    a=p.parse_args()
    print(json.dumps(run(a.archive,a.seed,a.coarse_support,a.checkpoints),indent=2,allow_nan=False))
