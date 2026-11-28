"""Moment-only bounds on high energy stored in collisionless cohort receipts."""

import argparse
import json
from pathlib import Path

import numpy as np

from frw_exterior_restart_reference import import_restart
from three_species_decay_onsager_reference import M, m, mp
from warm_light_occupation_stage_audit_reference import digest, read
from warm_prompt_reentry_shape_audit_reference import get_mom


def bounds(data,cut):
    a,d,j,s,bd,n,e=data
    r=np.hypot(d*np.sqrt((1-a)*(1+a)),a*e)
    floor=np.where(s==0,d,np.minimum(bd,e))
    ceiling=np.where(s==0,np.maximum(bd,e),np.inf)
    rest=d*np.sqrt((1-a)*(1+a))
    red=lambda x:np.hypot(rest,a*x)
    g=lambda x:np.maximum(red(x)-cut,0.)
    t=g(floor)
    finite=np.isfinite(ceiling)
    # For finite birth support the graph is below its chord. For a
    # semi-infinite birth interval use the derivative upper a.
    top=t+a*(e-floor)
    hi=ceiling[finite]
    lower=floor[finite]
    top[finite]=t[finite]+np.divide(g(np.where(finite,ceiling,floor))[finite]-t[finite],
                                    hi-lower,out=np.zeros_like(hi),where=hi>lower)*(e[finite]-lower)
    jensen=np.maximum(r-cut,0.)

    assert cut>np.max(d,initial=0.)
    x=np.sqrt((cut-rest)*(cut+rest))/a
    # The broken affine envelope is concave: slope on the left
    # is cut/(x-floor) >= a when floor<x. To the right the
    # redshifting energy increases at most a per unit birth energy.
    dominated=(floor>=x)
    left=np.divide(cut*(e-floor),x-floor,out=np.zeros_like(x),where=~dominated)
    right=cut+a*(e-x)
    tagged=np.where(dominated,red(floor)+a*(e-floor),np.where(e<=x,left,right))
    tagged=np.where(ceiling<x,0.,tagged)
    assert np.all(jensen<=top+5e-15*np.maximum(e,1e-280))
    assert np.all(tagged>=0.)
    return np.stack((n*jensen,n*top,n*tagged))


def independent():
    rng=np.random.default_rng(195)
    worst=0.
    for _ in range(1200):
        a=rng.uniform(.1,1.)
        d=rng.uniform(.001,.06)
        s=rng.integers(2)
        bd=d+rng.uniform(.001,.65)
        p=np.r_[rng.uniform(d,bd,18) if s==0 else bd+np.exp(rng.uniform(-9.,9.,18))]
        w=rng.random(len(p));w/=w.sum()
        x=np.array([a,d,0.,s,bd,1.,w@p])[:,None]
        cut=rng.uniform(.08,1.)
        B=bounds(x,cut)[:,0]
        E=np.hypot(d*np.sqrt((1.-a)*(1.+a)),a*p)
        obs=w@np.maximum(E-cut,0.)
        tail=w@(E*(E>=cut))
        worst=max(worst,float(max(B[0]-obs,obs-B[1],tail-B[2])))
    assert worst<2e-12
    return worst


def run(prior,paths):
    root=read(prior)["branches"][-1]
    raw=root["energy_upper_GeV4"]
    cut=np.hypot(M,1./root["entropy_expansion_from_original_anchor"])
    out=[]
    for path in paths:
        q=read(path)
        S=q["snapshot"]
        assert q["snapshot_sha256"]==digest(S)
        ex=import_restart(S["positive_exterior_restart"])
        P=np.asarray(S["lab_momentum_edges_GeV"])
        R=get_mom(P,[S["entropy_variable_"+b] for b in ("scalar","fermion")])[2].sum()
        moments=[]
        arr=ex.array()
        for C in (.10,.30,.50,cut):
            B=bounds(arr,C)
            assert np.max(B.sum(axis=1))<=ex.upper_bound(1.)*1.000000000001
            a=np.column_stack([B[:,arr[1]==d].sum(axis=1) for d in (mp,m)])
            moments.append(dict(threshold_GeV=C,absolute_GeV4=a.tolist(),
                                over_retained_energy=(a/R).tolist(),over_late_full_allowance=(a/raw).tolist()))
        out.append(dict(snapshot_sha256=q["snapshot_sha256"],receipt_count=arr.shape[1],
                        retained_energy_GeV4=R,cohorts_energy_lower_upper_GeV4=[ex.energy(1.),ex.upper_bound(1.)],
                        rows=moments))
    return dict(scope="Stored collisionless receipts only. The unknown birth-energy distributions are arbitrary "
                      "within each mean and support with subsequent ballistic redshift. No bound on old "
                      "unarchived energy, renewed growth or reacting re-entry.",
                columns=["excess_energy_lower","excess_energy_upper","energy_above_cut_upper"],
                species=["scalar","light"],synthetic_enclosure_max_failure=independent(),archives=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("prior",type=Path)
    p.add_argument("archives",nargs="+",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.prior,a.archives),indent=2,allow_nan=False))
