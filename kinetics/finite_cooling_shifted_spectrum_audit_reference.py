"""Same-metric comparison across the tiny displacement of the count roots."""
import argparse
import json
from pathlib import Path

import numpy as np

from finite_cooling_matched_history_los_reference import KEY, sha


def load(p):
    return json.loads(p.read_text())


def compare(prompt,mix,recoil,old,new):
    P,M,R=map(load,(prompt,mix,recoil))
    assert P["eos_sha256"]==M["eos_sha256"]==R["eos_sha256"]
    assert M["pairs_sha256"]==R["source_sha256"] and M["prompt_sha256"]==P["distribution_sha256"]
    rows=[]
    for label,A,B in zip(("instant","ratio5"),old,new,strict=True):
        a,b=load(A),load(B)
        assert a["z"]==b["z"]==20.
        assert a["eos_sha256"]==b["eos_sha256"]==M["eos_sha256"]
        pr=next(x["sha256"] for x in P["components"] if (x["case"],x["species"],x["sigma_ln_q"])==(label,"scalar",.14))
        mj=next(x["sha256"] for x in M["components"] if (x["case"],x["sigma_ln_q"])==(label,.14))
        re=R["components"][label]["table_sha256"]
        assert [mj,pr,re]==b["spectrum_sha256"]
        original={round(x["k_h_per_Mpc"],10):x for x in a["results"]}
        assert [round(x["k_h_per_Mpc"],10) for x in b["results"]]==[4.,10.,22.]
        data=[]
        for q in b["results"]:
            z=original[round(q["k_h_per_Mpc"],10)]
            for key in ("metric_sha256","original_native_metric_sha256","initial_k_tau","tau_final_Mpc"):
                assert z[key]==q[key]
            v0=np.array([x[0] for x in z["cases"][-1][KEY]])
            v1=np.array([x[0] for x in q["cases"][-1][KEY]])
            defect=float(abs(v1[0]-q["cases"][-1]["summed_prompt_recoil_with_full_IC"]))
            assert abs(v1[0]-v0[0])<.01 and defect<3e-7
            data.append(dict(k_h_Mpc=q["k_h_per_Mpc"],previous=v0.tolist(),shifted=v1.tolist(),change=(v1-v0).tolist(),
                             metric_sha256=q["original_native_metric_sha256"],additivity_defect=defect,
                             shifted_refinement=q["max_time_and_momentum_refinement_difference"]))
        rows.append(dict(history=label,probes=data,raw_hashes=[sha(A),sha(B)]))
    return dict(rows=rows,table_catalog_hashes=[sha(x) for x in (prompt,mix,recoil)],
        scope="Matched frozen scalar density tracer at the former and displaced count roots. Pair recoil before the same 1 MeV anchor, no post-anchor pair reshaping, stable distribution depletion, metric response, or total transfer.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for x in ("prompt","mix","recoil"):
        p.add_argument(x,type=Path)
    p.add_argument("--old",type=Path,nargs=2,required=True)
    p.add_argument("--new",type=Path,nargs=2,required=True)
    a=p.parse_args()
    print(json.dumps(compare(a.prompt,a.mix,a.recoil,a.old,a.new),indent=2,allow_nan=False))
