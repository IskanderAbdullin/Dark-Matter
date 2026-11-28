"""Check that saved positive tails contribute once across warm stage boundaries."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from coupled_warm_pair_continuation_reference import run
from frw_exterior_cohort_reference import Exterior
from frw_exterior_restart_reference import export_restart, import_restart
from scalar_pair_expanding_continuation_reference import moments
from three_species_decay_onsager_reference import m, mp


def read(p):
    return json.loads(p.read_text())


def snap(p):
    return read(p)["snapshots"][-1]


def ext(s):
    return import_restart(s["positive_exterior_restart"])


def public_count(s):
    # The conserved total includes each archived cohort exactly once.
    return np.sum(s["represented_species_count_a3"]) + np.sum(ext(s).array()[5])*s["scale_factor"]**3


def shell(s):
    return moments(np.asarray(s["lab_momentum_edges_GeV"]),
                   [np.asarray(s[k]) for k in ("entropy_variable_scalar","entropy_variable_fermion")])


def tracked(eos,s,stop,steps):
    r={}
    z=run(eos,s,steps,stop,checkpoint=r,track_exterior=True)
    for key in ("count_ledger","energy_ledger","collision_number_ledger"):
        assert abs(z[key])<2e-10
    return z,r


def direct_split(eos,s):
    assert abs(s["T_visible_GeV"]-.001)<1e-12 and "positive_exterior_restart" not in s
    x=Exterior()
    upper=s["lab_momentum_edges_GeV"][-1]
    p=2*upper
    n=.02*s["represented_species_count_a3"][0]/s["scale_factor"]**3
    E=np.hypot(p,mp)
    x.add(1.,mp,3,1,upper,n,n*p*p/(E+mp))
    s["positive_exterior_restart"]=export_restart(x)
    end=1.005
    half=np.exp(.5*np.log(end))
    A,a=tracked(eos,s,end,2)
    _,b=tracked(eos,s,half,1)
    B,c=tracked(eos,b,half,1)
    X,Y=shell(a),shell(c)
    diff=float(np.max(np.sum(abs(X-Y),axis=2)/np.sum(abs(X),axis=2)))
    counts=max(abs(public_count(z)/public_count(s)-1.) for z in (a,c))
    ex=[ext(z) for z in (a,c)]
    energy=abs(ex[0].energy(1.)/ex[1].energy(1.)-1.)
    hi=abs(ex[0].upper_bound(1.)/ex[1].upper_bound(1.)-1.)
    trace=abs(ex[0].trace/ex[1].trace-1.)
    dt=abs(a["T_visible_GeV"]/c["T_visible_GeV"]-1.)
    assert max(diff,counts,energy,hi,trace,dt)<2e-10
    assert A["exterior"]["predecessor_cohort_count"]==1 and B["exterior"]["predecessor_cohort_count"]>1
    return dict(direct_and_split_energy_ledger_defects=[A["energy_ledger"],B["energy_ledger"]],
                direct_split_T_relative=float(dt),direct_split_shell_l1_relative=diff,
                total_number_conservation_defect=float(counts),
                relative_differences_in_energy_upper_and_trace=[float(energy),float(hi),float(trace)],
                direct_exit={k:a[k] for k in ("T_visible_GeV","scale_factor","represented_species_count_a3")})


def producer_cross(producer, p_shape, base_producer, base_shape, next_run, next_shape):
    a,b=read(producer),read(base_producer)
    assert a==b and producer.read_bytes()==base_producer.read_bytes()
    s,t=snap(p_shape),snap(base_shape)
    assert t=={k:v for k,v in s.items() if k!="positive_exterior_restart"}
    q=a["rows"][-1]
    end=snap(next_shape)
    r=read(next_run)["rows"][-1]
    E=ext(s)
    N=E.array()[5]*s["scale_factor"]**3
    by_mass=np.array([N[E.array()[1]==mi].sum() for mi in (mp,m)])
    tails=np.asarray(q["first_pulse_count_and_total_energy_tails"])+np.asarray(q["comoving_source_channel_count_and_total_energy_tails"])
    v=tails[:,0].sum(axis=1)
    expected=np.array([v[2],v[0]+v[1]])+np.asarray(q["pair_exterior_counts_relative_to_final_retained"])*s["represented_species_count_a3"]
    E4=E.energy(1.)*s["scale_factor"]**4
    peer=q["visible_bath_feedback"]
    defects=np.r_[by_mass/expected-1.,E4/peer["exterior_a4_energy_one_mean_per_cohort"]-1.,
                  E.upper_bound(1.)*s["scale_factor"]**4/peer["exterior_a4_energy_distribution_independent_upper_bound"]-1.]
    new=ext(end)
    assert len(E.cohorts)==peer["exterior_cohorts"]==r["exterior"]["predecessor_cohort_count"]>0
    assert len(new.cohorts)>len(E.cohorts)
    assert r["exterior"]["initial_relative_a4_exterior_energy"]==E.energy(1.)
    joint=abs(public_count(end)/public_count(s)-1.)
    assert max(joint,np.max(abs(defects)),abs(r["energy_ledger"]))<2e-10
    return dict(producer_visible_T_GeV=s["T_visible_GeV"],continuation_visible_T_GeV=end["T_visible_GeV"],
                producer_count_of_positive_cohorts=len(E.cohorts),ending_count_of_positive_cohorts=len(new.cohorts),
                previous_a4_energy_EOS_units=E4,exterior_count_by_scalar_stable=by_mass.tolist(),
                comparison_to_producer_counts_energy_and_upper_bound_relative=defects.tolist(),
                joint_number_conservation_defect=joint,finished_energy_ledger=r["energy_ledger"])


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for key in ("eos","cold_shape","producer","producer_shape","legacy_producer","legacy_shape","bridge","bridge_shape"):
        p.add_argument(key,type=Path)
    a=p.parse_args()
    path=[a.eos,a.cold_shape,a.producer,a.producer_shape,a.legacy_producer,a.legacy_shape,a.bridge,a.bridge_shape]
    P=producer_cross(*path[2:])
    print(json.dumps(dict(producer_cross=P,synthetic_direct_split=direct_split(a.eos,snap(a.cold_shape)),
        input_sha256={v.name:hashlib.sha256(v.read_bytes()).hexdigest() for v in path},
        scope="Controlled carryover of known positive exterior cohorts in expansion only; no old-tail collision occupations or re-entry, and no inference about unarchived histories or terminal abundance."),indent=2,allow_nan=False))
