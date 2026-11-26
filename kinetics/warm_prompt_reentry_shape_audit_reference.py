"""Compare redistributed tails with the tail already present in q32."""

import argparse
import json
from pathlib import Path

import numpy as np

from joint_prompt_two_moment_evolution_reference import Collision
from joint_resonant_frw_split_reference import pack
from late_archived_shape_continuation_reference import Late,collide
from scalar_pair_expanding_continuation_reference import moments,MASS
from three_species_decay_onsager_reference import width
from warm_light_occupation_stage_audit_reference import digest,read


def pad(s,p):
    x=np.asarray(s["lab_momentum_edges_GeV"])
    assert len(p)>=len(x) and np.allclose(p[:len(x)],x,rtol=3e-14,atol=0.)
    ans=[]
    for key in ("entropy_variable_scalar","entropy_variable_fermion"):
        a=np.asarray(s[key])
        a=np.pad(a,((0,0),(0,len(p)-len(x))))
        a[0,len(x)-1:]=1000.
        ans.append(a)
    return ans


def get_mom(p,th):
    a=moments(p,th)
    counts=a[:,0]
    kin=a[:,1]
    return counts,kin,kin+MASS[:,None]*counts


def features(counts,K,E,core,old_wide):
    return dict(outside_q16_N_fraction=(counts[:,core:].sum(axis=1)/counts.sum(axis=1)).tolist(),
                outside_q16_energy_fraction=(E[:,core:].sum(axis=1)/E.sum(axis=1)).tolist(),
                outside_q32_N_fraction=(counts[:,old_wide:].sum(axis=1)/counts.sum(axis=1)).tolist(),
                outside_q32_energy_fraction=(E[:,old_wide:].sum(axis=1)/E.sum(axis=1)).tolist(),
                outside_q16_shell_number_weights=(counts[:,core:]/counts[:,core:].sum(axis=1)[:,None]).tolist(),
                outside_q16_shell_kinetic_weights=(K[:,core:]/K[:,core:].sum(axis=1)[:,None]).tolist())


def run(folder):
    D=folder
    F=[read(D/f"warm-light-occupation-next-q{q}-tape-2026-11-26.json.gz")[0] for q in (16,32)]
    A,B=[z["snapshot"] for z in F]
    assert all(z["snapshot_sha256"]==digest(z["snapshot"]) for z in F)
    core=len(A["lab_momentum_edges_GeV"])-1
    P=np.asarray(B["lab_momentum_edges_GeV"])
    old_wide=len(P)-1
    p=np.r_[P,P[-1]*np.geomspace(1.,2.,6)[1:]]
    H=read(D/"warm-near-rest-light-inverse-decay-2026-11-26.json.gz")["cases"][0]["stages"][0]["fine"]["H_over_rest_decay_width"]*width
    a=features(*get_mom(p,pad(B,p)),core,old_wide)
    rows=[]
    for tag,S in (("q16_initial",A),("q32_initial",B)):
        th=pad(S,p); v=pack(th)
        N0,_,E0=get_mom(p,th)
        seed=Late(Collision(p),H,[],True,local_jac=True)
        vals=[]
        times=[0.,1e-8,5e-8,2e-7,1e-6,5e-6]
        for x,y in zip(times[:-1],times[1:],strict=True):
            v,r,ok,no=collide(seed,v,y-x,S["T_visible_GeV"],continuation=True,balanced=True)
            N,K,E=get_mom(p,seed.arrays(v))
            er=[float(np.max(abs(N.sum(axis=1)/N0.sum(axis=1)-1.))),float(abs(E.sum()/E0.sum()-1.))]
            assert max(er)<2e-12
            vals.append(dict(elapsed_times_H=y,solver=[r,ok,no],invariant_defects=er,
                             **features(N,K,E,core,old_wide)))
        data=vals[-1]
        comp={k:(np.asarray(data[k])/np.asarray(a[k])-1).tolist() for k in
              ("outside_q16_N_fraction","outside_q16_energy_fraction")}
        shape={k:np.sum(abs(np.asarray(data[k])-np.asarray(a[k])),axis=1).tolist() for k in
               ("outside_q16_shell_number_weights","outside_q16_shell_kinetic_weights")}
        rows.append(dict(initial_snapshot=tag,rows=vals,
                         terminal_fraction_relative_changes_from_actual_q32_archive=comp,
                         terminal_shape_L1_from_actual_q32_archive=shape))
    return dict(scope="Fixed clock and physical momenta; isolated reacting prompt operator only. "
                      "Actual q32 checkpoint already contains a represented q16-to-q32 population. "
                      "Padding is vacuum outside each snapshot's own retained domain. Similar moments "
                      "or shapes after a short test do not establish a common history.",
                input_snapshot_sha256s=[z["snapshot_sha256"] for z in F],q64_edges_GeV=p.tolist(),
                actual_q32_before_run=a,controls=rows)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tapes",type=Path,default=Path("results"))
    print(json.dumps(run(parser.parse_args().tapes),indent=2,allow_nan=False))
