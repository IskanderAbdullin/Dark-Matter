"""Bound the high-momentum share of one-way exterior receipts from recorded moments."""

import argparse
import json
from pathlib import Path

import numpy as np

from frw_exterior_restart_reference import import_restart
from warm_light_occupation_stage_audit_reference import digest, read
from three_species_decay_onsager_reference import M, m


def tail(side, floor, n, mean, gate):
    """Sharp separate upper bounds on number and birth energy above a threshold."""
    if side == 0:
        assert floor<gate
        return 0.,0.
    if floor>=gate:
        return n,n*mean
    surplus=max(mean-floor,0.)
    return min(n,n*surplus/(gate-floor)), min(n*mean,n*gate*surplus/(gate-floor))


def test():
    # Two-point measure saturates both bounds simultaneously if its upper
    # point is at the gate and its lower point is at the known floor.
    f,g=2.,6.
    x,y=tail(1,f,4.,3.,g)
    assert x==1. and y==6.
    assert tail(0,2.,4.,1.,g)==(0.,0.)
    assert tail(1,7.,4.,8.,g)==(4.,32.)
    rng=np.random.default_rng(32842)
    for _ in range(60):
        energy=f+rng.lognormal(1.,2.,50)
        w=rng.random(50)
        a,b=tail(1,f,w.sum(),(energy@w)/w.sum(),g)
        high=energy>=g
        assert w[high].sum()<=a*(1+1e-15)
        assert w[high]@energy[high]<=b*(1+1e-15)
    return dict(two_point_number=x,two_point_energy=y)


def run(D):
    controls=test()
    data=read(D/"warm-heavy-high-band-finite-daughter-support-2026-11-26.json")
    rows=[]
    for r in data["cases"]:
        q=r["q"]; P=r["assumed_equal_comoving_daughter_cap_float_GeV"]; H=r["high_cut_GeV"]
        assert P<H/2 and m<M/2
        A=read(D/f"warm-early-q{q}-time32-checkpoint-2026-11-26.json.gz")
        B=read(D/f"warm-observed-suboperator-q{q}-tape-2026-11-26.json.gz")["blocks"][-1]["checkpoint"]
        S,T=A["snapshot"],B["snapshot"]
        assert A["snapshot_sha256"]==digest(S) and B["snapshot_sha256"]==digest(T)
        b=T["scale_factor"]/S["scale_factor"]
        before=import_restart(S["positive_exterior_restart"]).array()
        after=import_restart(T["positive_exterior_restart"]).array()
        n0=before.shape[1]
        after[[0,5]]*= np.asarray([b,b**3])[:,None]
        np.testing.assert_allclose(before,after[:,:n0],rtol=3e-14,atol=0.)
        summary=[]
        for stage,set_ in [("inherited_tracked",after[:,:n0]),("new_tracked",after[:,n0:])]:
            bi,mass,_,side,bound,N,E=set_
            assert np.all(mass<=m)
            # Energy at each cohort's exact birth expansion for q=H/2.
            threshold=np.hypot(mass,(H/2)/bi)
            estimates=np.asarray([tail(*z) for z in zip(side,bound,N,E,threshold,strict=True)])
            outer=side==1
            assert np.any(outer) and np.all(threshold[outer]>bound[outer])
            total=N[outer]@E[outer]
            high=estimates[:,1].sum()
            summary.append(dict(stage=stage,cohorts=len(N),outer_cohorts=int(outer.sum()),
                                outer_birth_energy_GeV4=total,
                                eligible_number_a3_GeV3=estimates[:,0].sum(),
                                eligible_birth_energy_a3_GeV4=high,
                                upper_share_of_tracked_outer_birth_energy=high/total,
                                largest_mean_birth_energy_to_threshold_ratio=float(np.max(E/threshold))))
        rows.append(dict(q=q,b_final=b,high_H_cut_GeV=H,
                         gross_ledger_input_tag_q_min_GeV=H/2,summary=summary))
    return dict(scope="Receipt-only moment bounds for stored exterior cohorts, without "
                      "resolving them as point masses. The resulting energy is a one-use "
                      "allowance for tracked cohorts only if momentum is not re-increased "
                      "after injection. No bound on actually reacting, repeatedly used, "
                      "untracked or newly accelerated inputs.",controls=controls,cases=rows)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tapes",type=Path,default=Path("results"))
    print(json.dumps(run(parser.parse_args().tapes),indent=2,allow_nan=False))
