"""Validate matched scripted observer tapes without treating trial states as accepted."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from warm_coupled_serial_restart_reference import bytes_of
from warm_light_occupation_stage_audit_reference import read
from warm_observed_suboperator_continuation_reference import metrics


def digest(z):
    return hashlib.sha256(bytes_of(z)).hexdigest()


def peaks(s):
    return metrics([s["entropy_variable_scalar"],s["entropy_variable_fermion"]])


def audit(tape,previous):
    history=read(tape)
    seed=read(previous)[-1]["snapshot"]
    parent_hash=digest({"snapshots":[seed]})
    base=seed["scale_factor"]
    rows=[]
    order=["redshift_1_input","redshift_1_output","pair_reconstruction_1",
          "pair_reconstruction_2","implicit_collision_output","pair_reconstruction_3",
          "pair_reconstruction_4","redshift_2_input","redshift_2_output"]
    for j,b in enumerate(history["blocks"],start=1):
        c,z=b["checkpoint"],b["observer"]
        assert digest(c)==z["checkpoint_sha256"]
        assert c["plan_sha256"]==z["plan_sha256"]==digest(c["plan"])
        assert c["plan"]["source_sha256"]==parent_hash
        assert c["blocks_finished"]==z["blocks_finished"]==j
        assert c["blocks_total"]==2
        a=z["samples"]; s=c["snapshot"]
        assert [r["operation"] for r in a]==order
        assert a[0]["entropy_array_sha256"]==peaks(seed)["entropy_array_sha256"]
        assert a[-1]["entropy_array_sha256"]==peaks(s)["entropy_array_sha256"]
        assert digest(s)==c["snapshot_sha256"]
        trial=np.asarray([r.get("trial_midpoint",False) for r in a])
        f=np.asarray([r["retained_light_peak"] for r in a])
        h=np.asarray([r["retained_scalar_peak"] for r in a])
        rows.append(dict(b_relative=s["scale_factor"]/base,T_GeV=s["T_visible_GeV"],
                         light_peak_start_and_end=[float(f[0]),float(f[-1])],
                         scalar_peak_start_and_end=[float(h[0]),float(h[-1])],
                         peak_all_samples=float(max(f)),peak_excluding_explicit_trials=float(max(f[~trial])),
                         midpoint_trials_count=int(sum(trial)),
                         retained_light_pair_quadratic_sign_margin_on_observed_nontrial_states=float(1-2*max(f[~trial])),
                         largest_mid_block_excess_over_both_endpoints=float(max(f)-max(f[0],f[-1])),
                         accepted=a[4]["microsteps_accepted"],rejected=a[4]["microsteps_rejected"]))
        seed=s
    ctrl=history.get("observer_free_first_checkpoint_sha256")
    if ctrl:
        assert ctrl==digest(history["blocks"][0]["checkpoint"])
    return dict(parent_source_sha256=parent_hash,hash_observer_free_first=ctrl,rows=rows)


def run(t16,t32,p16,p32):
    a,b=audit(t16,p16),audit(t32,p32)
    for x,y in zip(a["rows"],b["rows"],strict=True):
        assert abs(x["b_relative"]/y["b_relative"]-1)<1e-13
        assert abs(x["T_GeV"]/y["T_GeV"]-1)<1e-12
    return dict(scope="Instrumented explicit suboperator interfaces with fixed matched warm grids. "
                      "Two explicit pair midpoint evaluations marked trials. No sampling of all "
                      "internal nonlinear iterations, intervening times, heavy occupations, or "
                      "older exterior occupation shapes.",
                cases=[dict(q=16,**a),dict(q=32,**b)])


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("observer_q16","observer_q32","previous_q16","previous_q32"):
        p.add_argument(k)
    a=p.parse_args()
    print(json.dumps(run(*vars(a).values()),indent=2,allow_nan=False))
