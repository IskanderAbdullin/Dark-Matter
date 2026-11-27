"""Bound the post-archive direct-lepton vacuum source into the equality soft band."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from archived_collisionless_streaming_reference import Background
from decay_scalar_ir_relaxation_reference import limits
from prompt_scalar_cold_tail_envelope_reference import log_cap
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import M,P
from warm_late_cold_scalar_prompt_partition_reference import counts
from warm_light_occupation_stage_audit_reference import digest,read


def upper(bg,hi,lo,yield_scale,n):
    T=np.geomspace(hi,lo,n+1)
    scale=np.array([bg.scale(t) for t in T])
    p=.1*bg.T0/scale
    assert np.all(np.diff(scale)>0.) and np.all(np.diff(p)<0.) and np.all(p<P)
    minimum_parent=np.array([limits(q)[1] for q in p])
    assert np.all(np.diff(minimum_parent/T)>0)
    # Every bin uses the early upper source and lower photon-only H*s
    # at the late endpoint. This upper remains one-sided for an EOS
    # through electron annihilation and avoids integrating a fitted rate.
    src=[log_cap(t,q)[0] for t,q in zip(T[:-1],p[:-1],strict=True)]
    photon_H=np.sqrt(np.pi**2/45)*T[1:]**2/MPL
    photon_s=4*np.pi**2/45*T[1:]**3
    terms=np.log(2*np.diff(np.log(scale)))+np.asarray(src)-np.log(photon_H*photon_s)
    z=float(logsumexp(terms))
    return dict(temperature_panels=n,log10_vacuum_soft_yield_upper=z/np.log(10.),
                log10_upper_over_archive_cold_yield=(z-np.log(yield_scale))/np.log(10.))


def clock_monotonicity(bg):
    # Certify d(3 log(T)+log(gs))/dlog(T)>0 on every cubic piece,
    # not only at the integration grid used to form the envelope.
    a,b,c,_=np.asarray(bg.gs.c)
    h=np.diff(bg.gs.x)
    stationary=np.divide(-b,3*a,where=a!=0,out=np.zeros_like(a))
    y=np.clip(stationary,0,h)
    def val(t):
        return 3+3*a*t*t+2*b*t+c
    lo=float(min(np.min(val(np.zeros_like(h))),np.min(val(h)),np.min(val(y))))
    assert lo>0.
    return lo


def run(table,source):
    bg=Background(table)
    clock=clock_monotonicity(bg)
    s=read(source)
    S=s["snapshot"]
    assert s["snapshot_sha256"]==digest(S)
    start=S["T_visible_GeV"]
    ai=bg.scale(start)
    p=.1*bg.T0/ai
    n=counts(S,np.array([p]),24)[0]
    entropy=2*np.pi**2/45*np.exp(bg.gs(np.log(start)))*start**3
    Y=float(n/entropy)
    lines=[]
    for T0 in (start,.010442,.012,.025,.05):
        a=bg.scale(T0)
        cut=.1*bg.T0/a
        scale=limits(cut)[1]
        later=.0001
        # Use dyadic pairs to enforce nested partition monotonicity.
        estimates=[upper(bg,T0,later,Y,n) for n in (32,128,512)]
        assert all(x["log10_vacuum_soft_yield_upper"]>=y["log10_vacuum_soft_yield_upper"] for x,y in zip(estimates,estimates[1:]))
        lines.append(dict(T_max_GeV=T0,T_min_GeV=later,p_over_T_at_start=cut/T0,
                            earliest_one_parent_plus_partner_at_rest_over_T=float((scale+M)/T0),
                            refinements=estimates))
    return dict(scope="Direct thermal charged-lepton vacuum source only, with immediate two-body decay, "
                      "no initially present heavy, hadronic production, inverse reactions or later scattering. "
                      "For each row integrate from that start to 0.1 MeV. The denominator is only the "
                      "same known archive's held equality-cold yield, not its earlier-temperature history.",
                archive=s["snapshot_sha256"],minimum_log_entropy_slope_on_EOS_cubics=clock,comparison_archive_held_equality_cold_yield=Y,rows=lines)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("archive",type=Path)
    args=p.parse_args()
    print(json.dumps(run(args.eos,args.archive),indent=2,allow_nan=False))
