"""Audit matched aligned closed prompt time ladders without inferring an unsampled turning time."""

import argparse
import json
from pathlib import Path

import numpy as np

from warm_light_occupation_stage_audit_reference import digest,read


CLOCK=np.r_[0.,np.geomspace(1e-11,1e-6,6)]


def curve(D,n,count):
    rows=[]
    src=None
    shape=None
    for i,t in enumerate(CLOCK[:count+1]):
        v=read(D/f"step-{i}.json")
        assert v["time"]==t and v["hash"]==digest(v["snapshot"])
        plan=v["plan"]
        assert plan["substeps"]==n and plan["artificial_empty_shell_filling"]==0. and plan["positive_predictor"]
        if src is None:
            src=plan
            shape=v["snapshot"]
            assert src["seed_sha256"]==v["hash"]
        assert plan==src
        assert v["snapshot"]["lab_momentum_edges_GeV"]==shape["lab_momentum_edges_GeV"]
        assert v["snapshot"]["T_visible_GeV"]==shape["T_visible_GeV"]
        assert np.max(abs(np.asarray(v["moment_defects"])))<6e-12
        W=np.asarray(v["excess_over_initial_energy"])
        Q=np.asarray(v["independent_closed_prompt_gross_and_signed"])
        assert len(W)==2 and np.min(W)>=0.
        assert abs(Q[3]-(Q[5]-Q[4]))<5e-14*max(Q[1],Q[2],1e-230)
        assert v["cut_moment_error_over_rho"]<5e-11*max(W.sum(),1e-60)
        rows.append(dict(state_sha256=v["hash"],time=float(t),excess_energy=W.tolist(),
                         signed_native_flux=v["native_signed_cut_flux"],angular_check=Q[3],
                         signed_flux_discrepancy_over_gross=float(abs(v["native_signed_cut_flux"]-Q[3])/
                                                          max(Q[1]+Q[2],1e-240))))
    sign=np.sign([r["signed_native_flux"] for r in rows])
    edges=[[rows[i]["time"],rows[i+1]["time"]] for i in range(len(sign)-1) if sign[i]*sign[i+1]<0]
    return src,rows,edges


def run(coarse,fine,count):
    a,A,b=curve(coarse,3,count)
    c,C,d=curve(fine,6,count)
    assert a["hard_GeV"]==c["hard_GeV"] and a["seed_sha256"]==c["seed_sha256"]
    assert a["eos_sha256"]==c["eos_sha256"]
    w=np.array([[sum(f["excess_energy"]) for f in Q] for Q in (A,C)])
    ref=float(w.max())
    assert ref>0.
    x,t=w[:,0],w[:,-1]
    return dict(scope="Identical exact-vacuum initial moments and cut-aligned finite pad; isolated prompt at fixed "
                      "clock. Positive and negative samples bracket at least one sign change, not the first/unique turning time.",
                hard_packet_GeV=a["hard_GeV"],completed_intervals=count,
                two_time_ladders=[A,C],opposite_sign_brackets=[b,d],
                max_excess_difference_over_largest_sampled_excess=float(np.max(abs(w[0]-w[1]))/ref),
                terminal_minus_initial_excess_over_initial_dark_energy=(t-x).tolist())


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("coarse",type=Path)
    p.add_argument("fine",type=Path,nargs="?")
    p.add_argument("--intervals",type=int,default=6)
    a=p.parse_args()
    assert 0<=a.intervals<len(CLOCK)
    if a.fine is not None:
        result=run(a.coarse,a.fine,a.intervals)
    else:
        src,rows,sign=curve(a.coarse,3,a.intervals)
        result=dict(scope="One fixed-clock cut-aligned prompt-only ladder, no time refinement.",
                    hard_packet_GeV=src["hard_GeV"],completed_intervals=a.intervals,
                    opposite_sign_brackets=sign,rows=rows)
    print(json.dumps(result,indent=2,allow_nan=False))
