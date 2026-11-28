"""Isolate the mixed time/support difference at matched corrected stages."""

import argparse
import json
from pathlib import Path

import numpy as np

from warm_history_lower_matched_audit_reference import measures
from warm_light_occupation_stage_audit_reference import digest, read
from warm_multishell_corrected_support_audit_reference import audit


def run(eos, seed32, narrow24, narrow48, wide24, wide48, count):
    A=audit(eos,seed32,narrow24,wide24,count)
    B=audit(eos,seed32,narrow48,wide48,2*count)
    root=read(seed32)[-1]
    s0=root["snapshot"]
    X=measures(s0)*s0["scale_factor"]**3
    n=X.shape[-1]
    norm=X.sum(axis=-1)
    init=[root,read(wide24/"guide-start.json"),root,read(wide48/"guide-start.json")]
    assert init[1]["snapshot_sha256"]==init[3]["snapshot_sha256"]
    frames=[]
    for i in range(count+1):
        W=(A["frames"][i],B["frames"][2*i])
        q=np.array([v["T_GeV"] for f in W for v in f["branches"]])
        np.testing.assert_allclose(q,q[0],rtol=2e-13)
        objects=[read(D/f"step-{j}.json") for D,j in
                 ((narrow24,i),(wide24/"q64",i),(narrow48,2*i),(wide48/"q64",2*i))] if i else init
        res=[]
        p0=None
        for k,R in enumerate(objects):
            s=R["snapshot"]
            assert digest(s)==R["snapshot_sha256"]==W[k//2]["branches"][k%2]["snapshot_sha256"]
            p=np.asarray(s["lab_momentum_edges_GeV"])[:n+1]
            if p0 is None:
                p0=p
            np.testing.assert_allclose(p,p0,rtol=1.3e-13,atol=0.)
            res.append(measures(s)[:,:,:n]*s["scale_factor"]**3)
        a,b,c,d=res
        nt,wt=c-a,d-b
        ds0,ds1=b-a,d-c
        I=wt-nt
        assert np.all(abs(I-(ds1-ds0))<=1e-10*abs(ds1-ds0)+norm[:,:,None]*1e-15)
        cold=[x["first_cell_occupation"] for t in W for x in t["branches"]]
        def l1(x):
            return (abs(x).sum(axis=-1)/norm).tolist()
        frames.append(dict(coarse_index=i,T_GeV=float(q[0]),time_L1_q32=l1(nt),time_L1_q64=l1(wt),
                           support_L1_24=l1(ds0),support_L1_48=l1(ds1),
                           signed_time_support_interaction_L1=l1(I),cold_average=cold,
                           signed_cold_interaction=(cold[3]-cold[1])-(cold[2]-cold[0])))
    return dict(scope="Same-clock four-run difference of differences on the shared cells. "
                      "No bound on omitted momentum space or time extrapolation outside tested frames.",
                metric_order=["number","full_energy","pressure"],
                normalization="initial q32 represented per-species moment",
                branch_order=["q32-24","q64-24","q32-48","q64-48"],
                shared_coarse_intervals=count,frames=frames)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for x in ("eos","seed32","narrow24","narrow48","wide24","wide48"):
        p.add_argument(x,type=Path)
    p.add_argument("--steps",type=int,default=24)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.seed32,a.narrow24,a.narrow48,a.wide24,a.wide48,a.steps),indent=2,allow_nan=False))
