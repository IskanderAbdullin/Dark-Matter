"""Resolve reversible prompt crossings, not just the net exterior fraction."""

import argparse
import json
from pathlib import Path

import numpy as np

from joint_prompt_two_moment_evolution_reference import Collision
from joint_resonant_frw_split_reference import pack
from late_archived_shape_continuation_reference import Late,collide
from three_species_decay_onsager_reference import width
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom,pad


def modes(model,th,ncut):
    c=model.k
    f=c.c[0].f(th[0],c.ph,c.I)
    g=c.c[1].f(th[1],c.ch,c.J)
    A=c.pre*f*g
    B=c.pre*(1+f)*(1-g)
    def r(q,mask=None):
        return np.bincount(c.y,weights=q if mask is None else q*mask,minlength=c.n)
    den=r(A+B)
    factor=np.divide(c.w,den,out=np.zeros_like(den),where=den>0.)
    net=c.rate(th)[:,0]
    output=[]
    for j,label in enumerate((c.I,c.J)):
        ext=label>=ncut
        a,b=r(A,ext),r(B,ext)
        inward=float(factor@(a*r(B,~ext)))
        outward=float(factor@(b*r(A,~ext)))
        # Check against the already verified two-moment projection; compare
        # the difference on the gross scale, which stays finite at balance.
        direct=float(net[j,ncut:].sum())
        assert abs((outward-inward)-direct)<3e-12*max(outward+inward,1e-280)
        output.append([outward,inward,direct])
    return np.asarray(output)


def run(D):
    ref=read(D/"warm-near-rest-light-inverse-decay-2026-11-26.json.gz")
    C=[read(D/f"warm-light-occupation-next-q{q}-tape-2026-11-26.json.gz")[0] for q in (16,32)]
    A,B=[k["snapshot"] for k in C]
    assert all(k["snapshot_sha256"]==digest(k["snapshot"]) for k in C)
    Pc,Pr=[np.asarray(s["lab_momentum_edges_GeV"]) for s in (A,B)]
    grids=np.r_[Pr,Pr[-1]*np.geomspace(1.,2.,6)[1:]]
    ncut=len(Pc)-1
    H=ref["cases"][0]["stages"][0]["fine"]["H_over_rest_decay_width"]*width

    def walk(s,mesh,extend,steps,subdivision=1):
        th=pad(s,mesh); v=pack(th)
        n_init=get_mom(mesh,th)[0]
        N0=n_init.sum(axis=1)
        initial_out=n_init[:,extend:].sum(axis=1)/N0
        model=Late(Collision(mesh),H,[],True,local_jac=True)
        rows=[]
        gross=np.zeros((2,2))
        def evaluate(V):
            a=model.arrays(V)
            N=get_mom(mesh,a)[0]
            q=modes(model,a,extend)/(H*N0[:,None])
            return N[:,extend:].sum(axis=1)/N0,q
        def tab(t,V):
            population,q=evaluate(V)
            clos=np.max(abs(gross[:,0]-gross[:,1]-(population-initial_out)))
            assert clos<2e-12
            rows.append(dict(elapsed_times_H=t,number_outside_edge_fraction=population.tolist(),
                             gross_out_in_and_direct_net_per_initial_number_per_Hubble=q.tolist(),
                             discrete_cumulative_gross_out_in_per_initial_number=gross.tolist(),
                             maximum_discrete_count_closure_error=float(clos)))
        tab(0.,v)
        for t,u in zip(steps[:-1],steps[1:],strict=True):
            h=(u-t)/subdivision
            for _ in range(subdivision):
                v,_,accepted,_=collide(model,v,h,s["T_visible_GeV"],continuation=True,balanced=True)
                # Track the right-end flux of each actual accepted backward
                # Euler step. The approximation is not a continuous-time integral.
                assert accepted==1
                gross+=h*evaluate(v)[1][:,:2]
            tab(u,v)
        return rows
    t=[0.,1e-8,5e-8,2e-7,1e-6,5e-6]
    first=walk(A,grids,ncut,t)
    second=walk(B,grids,ncut,[0.])
    third=walk(B,grids,len(Pr)-1,t)
    refined=walk(A,grids,ncut,t,2)
    finer=walk(A,grids,ncut,t,4)
    def integrated(rows):
        return np.asarray(rows[-1]["discrete_cumulative_gross_out_in_per_initial_number"])
    change=np.max(abs(integrated(refined)/integrated(finer)-1))
    assert change<.004
    assert min(np.array(first[0]["gross_out_in_and_direct_net_per_initial_number_per_Hubble"])[:,0])>1000.
    assert max(np.array(first[0]["gross_out_in_and_direct_net_per_initial_number_per_Hubble"])[:,1])==0.
    assert min(np.array(first[-1]["gross_out_in_and_direct_net_per_initial_number_per_Hubble"])[:,1])>0.
    return dict(scope="Instantaneous factorized effective prompt crossings. Each virtual parent "
                      "energy has inverse weights A and decay weights B, with outgoing-sector "
                      "crossing requiring inner inverse and outer decay and vice versa. "
                      "The backward-Euler ledger counts effective signed transitions consistently "
                      "with each actual step, not an exact continuous-time or labelled-particle history.",
                input_snapshot_sha256s=[k["snapshot_sha256"] for k in C],
                from_q16_vacuum_across_q16_edge=first,from_real_q32_across_q16_edge=second,
                from_real_q32_across_q32_edge=third,
                q16_initial_refined_2=refined,q16_initial_refined_4=finer,
                maximum_terminal_gross_relative_change_for_last_two_orders=float(change))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tapes",type=Path,default=Path("results"))
    print(json.dumps(run(p.parse_args().tapes),indent=2,allow_nan=False))
