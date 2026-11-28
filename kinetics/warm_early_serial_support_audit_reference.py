"""Independent q16/q32 inclusive count audit on the same short warm grid."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from three_species_decay_onsager_reference import m, mp
from warm_exterior_restart_audit_reference import ext
from warm_coupled_serial_restart_reference import bytes_of


def read(p):
    return json.loads(p.read_text())


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def rows(S,B,arch,source_sha):
    assert "positive_exterior_restart" not in S
    assert B["plan"]["source_sha256"]==source_sha
    assert B["plan_sha256"]==hashlib.sha256(bytes_of(B["plan"])).hexdigest()
    assert B["snapshot_sha256"]==hashlib.sha256(bytes_of(B["snapshot"])).hexdigest()
    assert B["blocks_finished"]==len(B["ledger"])
    v=np.asarray(S["represented_species_count_a3"]).copy()
    e=np.zeros(2)
    R=np.asarray(arch["retained_fixed_Omega_h2"])
    X=np.asarray(arch["exterior_fixed_Omega_h2"])
    F=R/v
    assert abs(F[0]/mp/(F[1]/m)-1.)<1e-12
    u=v.copy()
    out=[]
    for l in B["ledger"]:
        a=np.asarray(l["endpoint"]["retained_counts_relative_to_start"])
        b=np.asarray(l["endpoint"]["inclusive_counts_relative_to_start"])
        diff=(b-a)*v
        # Recorded per-block relative count sums lose signed rounding
        # below ~1e-16 of the retained density. Check against the
        # independent stored positive-tail archive at the final state.
        assert np.all(diff>-2e-13*u)
        e+=diff
        v*=a
        assert abs((v+e).sum()/u.sum()-1.)<1e-9
        out.append(dict(b=l["relative_expansion"],T_GeV=l["T_visible_GeV"],
                        retained_fixed_count_omega=(F*v).tolist(),inclusive_fixed_count_omega=(F*(v+e)+X).tolist(),
                        retained_fraction_from_seed=(v/u).tolist(),new_tail_count_relative_to_initial=(e/u).tolist(),
                        energy_ledger=l["energy_ledger"],microsteps=[l["accepted_microsteps"],l["rejected_microsteps"]]))
    E=ext(B["snapshot"])
    a=B["snapshot"]["scale_factor"]
    q=E.array()
    direct=np.array([q[5,q[1]==mi].sum() for mi in (mp,m)])*a**3
    assert max(abs(direct-e)/u)<2e-13
    assert max(abs(v/B["snapshot"]["represented_species_count_a3"]-1.))<2e-10
    return out


def run(archive,seed16,cp16,seed32,cp32,allow_prefix=False):
    A=read(archive)["cases"]
    s,t=[read(p)["snapshots"][-1] for p in (seed16,seed32)]
    a,b=map(read,(cp16,cp32))
    assert abs(t["T_visible_GeV"]/s["T_visible_GeV"]-1.)<1e-12
    assert abs(t["scale_factor"]/s["scale_factor"]-1.)<1e-12
    q,r=[p["lab_momentum_edges_GeV"] for p in (s,t)]
    assert len(q)<len(r) and np.max(abs(np.array(q)-r[:len(q)]))<1e-12
    for p in ("schema","eos_sha256","outer_log_scale_grid","local_prompt_jacobian","reference_sha256"):
        assert a["plan"][p]==b["plan"][p]
    assert a["blocks_total"]==b["blocks_total"]==4
    if not allow_prefix:
        assert a["blocks_finished"]==b["blocks_finished"]==4
    P=rows(s,a,A[0],sha(seed16));Q=rows(t,b,A[1],sha(seed32))
    o=[]
    for p,q in zip(P,Q):
        assert p["b"]==q["b"] and abs(p["T_GeV"]/q["T_GeV"]-1)<1e-10
        l=np.asarray(q["inclusive_fixed_count_omega"])/p["inclusive_fixed_count_omega"]-1
        o.append(dict(b=p["b"],T_GeV=p["T_GeV"],q16=p,q32=q,relative_q32_minus_q16_inclusive=l.tolist()))
    return dict(complete=len(o)==4,matched=o,archive_sha256=[sha(p) for p in (archive,seed16,cp16,seed32,cp32)],
                scope="Independent short 12-MeV warm support control. Count old positive external tails once from the producer summaries; their historical distributions and gravity are absent in both seeds. Carry only newly resolved external cohorts in the warm clock. Neither a time-grid convergence result nor a late relic abundance.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for x in ("archive","seed16","cp16","seed32","cp32"):
        p.add_argument(x,type=Path)
    p.add_argument("--allow-prefix",action="store_true")
    a=p.parse_args()
    print(json.dumps(run(a.archive,a.seed16,a.cp16,a.seed32,a.cp32,a.allow_prefix),indent=2,allow_nan=False))
