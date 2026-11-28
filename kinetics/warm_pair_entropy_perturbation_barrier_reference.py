"""Frozen pair sign robustness under bounded errors in entropy variables."""

import argparse
import json
from pathlib import Path

import numpy as np

from warm_light_occupation_stage_audit_reference import digest, read


def bound(row, scalar_margin, light_margin):
    a,b=scalar_margin,light_margin
    g=row.get("entropy_gap_inside",row.get("represented_gap"))
    s=row["exterior_scalar_minimum_entropy"]
    gap=g-2*a-b
    if gap<=0 or (s is not None and s<=a):
        return dict(successful=False,slack=float(gap))
    loss=row.get("interior_loss_lower_GeV",row.get("negative_lower_GeV"))*np.exp(-b)*(-np.expm1(-gap))/(-np.expm1(-g))
    pos=row.get("exterior_gain_upper_GeV",row.get("positive_upper_GeV"))
    if s is not None:
        pos*=np.exp(2*a)*(np.expm1(-s)/np.expm1(-(s-a)))**2
    return dict(successful=bool(pos<loss),slack=float(gap),upper_over_lower=float(pos/loss))


def run(tracks):
    out={}
    for key,tape in tracks.items():
        out[key]={}
        for group, rows in (("frozen",tape["frozen"]["frames"]),("interpolant",tape["interpolant"]["brackets"])):
            tests=[]
            for a,b in [(1.,1.),(1.75,1.75),(2.,2.),(3.,3.),(3.5,3.5),(2.,6.)]:
                res=[bound(r,a,b) for r in rows]
                vals=[r.get("upper_over_lower") for r in res]
                tests.append(dict(scalar_entropy_sup_error=a,light_entropy_sup_error=b,
                                  all_boxes_certified=all(r["successful"] for r in res),
                                  minimum_inside_slack=min(r["slack"] for r in res),
                                  maximum_exterior_rate_ratio=None if any(v is None for v in vals) else max(vals)))
            out[key][group]=dict(count=len(rows),tests=tests)
    return dict(scope="Only declared sup entropy errors, fixed supports and vacuum incoming exterior. Actual errors and cap are not inferred.",cases=out)


def combine(frozen_path, i32, i64):
    tracks=read(frozen_path)
    for label,p in (("q32",i32),("q64",i64)):
        rows=read(p)
        assert all(x["successful"] for x in rows["brackets"])
        assert rows["brackets"][-1]["snapshot_sha256"][-1] == tracks[label]["frozen"]["frames"][-1]["snapshot_sha256"]
        tracks[label]["interpolant"]=rows
    return run(tracks)


def compare(a,b):
    entries=[read(v) for v in (a,b)]
    items=[x["snapshot"] for x in entries]
    assert all(digest(q)==x["snapshot_sha256"] for q,x in zip(items,entries,strict=True))
    A,B=items
    np.testing.assert_allclose(A["lab_momentum_edges_GeV"],B["lab_momentum_edges_GeV"],rtol=1e-13,atol=0.)
    assert abs(A["T_visible_GeV"]/B["T_visible_GeV"]-1.)<1e-13
    out=[]
    for k in ("scalar","fermion"):
        u=np.asarray(A["entropy_variable_"+k])-np.asarray(B["entropy_variable_"+k])
        out.append(float(np.max(abs(np.stack((u[0],u[0]+u[1]))))))
    return dict(snapshot_sha256=[x["snapshot_sha256"] for x in entries], T_GeV=A["T_visible_GeV"],
                maximum_scalar_light_difference=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("tracks",type=Path)
    p.add_argument("--interpolation",nargs=2,type=Path)
    p.add_argument("--compare",action="append",nargs=2,type=Path,default=[])
    a=p.parse_args()
    answer=run(read(a.tracks)) if a.interpolation is None else combine(a.tracks,*a.interpolation)
    answer["observed_accepted_numerical_discrepancies"]=[compare(*q) for q in a.compare]
    print(json.dumps(answer,indent=2,allow_nan=False))
