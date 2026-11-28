"""Compare one coarse q32 outer block with two half-size observed blocks."""
import argparse

import numpy as np

from warm_light_occupation_stage_audit_reference import read,digest
from frw_exterior_restart_reference import import_restart
from scalar_pair_expanding_continuation_reference import moments


def one(entry):
    C=entry["checkpoint"]; O=entry["observer"]
    assert digest(C)==O["checkpoint_sha256"]
    assert digest(C["snapshot"])==C["snapshot_sha256"]
    assert digest(C["plan"])==O["plan_sha256"]==C["plan_sha256"]
    A=O["samples"]
    assert len(A)==9 and A[3]["operation"]=="pair_reconstruction_2"
    assert not A[3]["trial_midpoint"]
    f=lambda x:x["retained_light_peak"]
    return dict(start=f(A[0]),end=f(A[-1]),first_kick=f(A[3]),peak=max(f(r) for r in A),
                max_excess=max(f(r) for r in A)-max(f(A[0]),f(A[-1])),
                first_pair_gain=f(A[3])-f(A[0]))


def receipt(initial,final):
    base=import_restart(initial["positive_exterior_restart"])
    end=import_restart(final["positive_exterior_restart"])
    b=final["scale_factor"]/initial["scale_factor"]
    a=base.array()
    c=end.array()
    n=a.shape[1]
    np.testing.assert_allclose(c[[0,5],:n]*(np.array([b,b**3]))[:,None],
                               a[[0,5]],rtol=3e-14,atol=0.)
    assert np.array_equal(c[[1,2,3,4,6],:n],a[[1,2,3,4,6]])
    z=c[:,n:]
    E=z[5]*z[6]*b**3
    assert abs(end.injected*b**4-base.injected-E@(z[0]*b))<1e-12*end.injected*b**4
    return float(sum(E[z[3]==1]))


def run(coarse,fine,previous):
    c=read(coarse)["blocks"][0]
    f=read(fine)["blocks"]
    assert len(f)==2 and all(z["checkpoint"]["plan"]["source_sha256"]==
                              c["checkpoint"]["plan"]["source_sha256"] for z in f)
    assert f[0]["observer"]["samples"][0]==c["observer"]["samples"][0]
    X,Y=c["checkpoint"]["snapshot"],f[-1]["checkpoint"]["snapshot"]
    S=read(previous)[-1]["snapshot"]
    assert digest({"snapshots":[S]})==c["checkpoint"]["plan"]["source_sha256"]
    ext=[receipt(S,t) for t in (X,Y)]
    assert abs(X["scale_factor"]/Y["scale_factor"]-1.)<1e-13
    assert abs(X["T_visible_GeV"]/Y["T_visible_GeV"]-1.)<1e-12
    p,q=[np.asarray(w["lab_momentum_edges_GeV"]) for w in (X,Y)]
    assert max(abs(p-q))<1e-15
    V=[]
    for s in (X,Y):
        V.append(moments(p,[np.asarray(s[k]) for k in ("entropy_variable_scalar","entropy_variable_fermion")]))
    A,B=V
    differences=(np.sum(abs(B-A),axis=2)/np.sum(A,axis=2)).tolist()
    q=one(c); r=[one(z) for z in f]
    # Use the same original block entrance for the total trajectory maximum.
    fine_peak=max(v["peak"] for v in r)
    fine_excess=fine_peak-max(q["start"],r[-1]["end"])
    assert fine_excess>=0. and q["max_excess"]>0. and fine_excess<q["max_excess"]
    return dict(scope="Same inherited q32 input and same combined endpoint: "
                      "coarse versus two half-size splitting blocks. Samples only "
                      "the declared explicit interface states, not the continuous solution.",
                endpoint_T_GeV=Y["T_visible_GeV"],coarse=q,fine=r,
                refined_combined_peak=fine_peak,
                refined_combined_excess_over_shared_interval_endpoints=fine_excess,
                fine_over_coarse_within_interval_excess=fine_excess/q["max_excess"],
                scalar_stable_shell_L1_number_and_kinetic=differences,
                gross_new_high_exterior_creation_coarse_fine_GeV4=ext,
                gross_new_high_exterior_creation_relative_change=ext[1]/ext[0]-1.,
                relative_total_counts=(np.sum(B[:,0],axis=1)/np.sum(A[:,0],axis=1)-1.).tolist())


if __name__=="__main__":
    import json
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("coarse")
    p.add_argument("fine")
    p.add_argument("previous")
    a=p.parse_args()
    print(json.dumps(run(a.coarse,a.fine,a.previous),indent=2,allow_nan=False))
