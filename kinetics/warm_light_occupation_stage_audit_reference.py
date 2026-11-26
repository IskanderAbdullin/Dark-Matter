"""Audit a new continued warm stage without changing the integrator."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.special import expit,roots_legendre

from frw_exterior_restart_reference import import_restart
from joint_visible_heavy_adiabatic_reference import parent_branch,snapshot_functions
from three_species_decay_onsager_reference import M
from warm_coupled_serial_restart_reference import bytes_of,sha
from warm_isolated_recycling_energy_envelope_reference import moments


def read(p):
    x=Path(p).read_bytes()
    return json.loads(gzip.decompress(x) if str(p).endswith(".gz") else x)


def digest(x):
    return hashlib.sha256(bytes_of(x)).hexdigest()


def stage_metrics(S):
    u,v=np.asarray(S["entropy_variable_fermion"])
    x,y=np.asarray(S["entropy_variable_scalar"])
    N,R=moments(S).sum(axis=0)
    F=float(expit(-min(np.min(u),np.min(u+v))))
    B=float(1/np.expm1(min(np.min(x),np.min(x+y))))
    ext=import_restart(S["positive_exterior_restart"])
    totals=float(N+np.sum(ext.array()[5]))*S["scale_factor"]**3
    f,ef,g,eg=snapshot_functions(S)
    cap=np.sqrt(max(0.,(ef[-1]+eg[-1])**2-M*M))
    trials=[]
    for order in (5,11):
        val=0.
        if cap>0:
            t=np.linspace(0,cap,60)
            a,b=(t[:-1,None]+t[1:,None])/2,np.diff(t)[:,None]/2
            x,w=roots_legendre(order)
            p=(a+b*x).ravel()
            e=np.hypot(M,p)
            H=parent_branch(e,f,g,ef,eg,angle=15)
            val=(2/np.pi**2)*np.sum((b*w).ravel()*p**2*e*H[:,0])/R
        trials.append(float(val))
    assert abs(trials[0]-trials[1])<4e-6*max(trials[-1],1e-19)
    return dict(b_absolute=S["scale_factor"],T_GeV=S["T_visible_GeV"],
                retained_light_peak=F,retained_scalar_peak=B,
                frozen_heavy_stationary_energy_over_retained=trials,
                total_tracked_comoving_number_GeV3=totals)


def verify(parent,tape):
    A=read(parent)
    B=read(tape)
    assert A["blocks_finished"]==A["blocks_total"]==32 and isinstance(B,list)
    S=A["snapshot"]
    assert A["snapshot_sha256"]==digest(S)
    target=digest({"snapshots":[S]})
    v=[stage_metrics(S)]
    for i,c in enumerate(B,1):
        assert c["blocks_finished"]==len(c["ledger"])==i
        assert c["blocks_total"]==4
        assert c["plan"]["source_sha256"]==target
        assert c["snapshot_sha256"]==digest(c["snapshot"]) and c["plan_sha256"]==digest(c["plan"])
        assert c["plan"]==B[0]["plan"]
        v.append(stage_metrics(c["snapshot"]))
    for w in v:
        w["b_relative_to_source"]=w["b_absolute"]/v[0]["b_absolute"]
    defect=max(abs(t["total_tracked_comoving_number_GeV3"]/v[0]["total_tracked_comoving_number_GeV3"]-1.) for t in v)
    assert defect<2e-9
    return dict(parent_sha256=sha(parent),tape_sha256=sha(tape),source_expected_sha256=target,
                tracked_number_relative_defect=defect,stages=v)


def run(p16,p32,t16,t32):
    a,b=verify(p16,t16),verify(p32,t32)
    assert len(a["stages"])==len(b["stages"])
    for x,y in zip(a["stages"],b["stages"],strict=True):
        assert abs(x["b_relative_to_source"]/y["b_relative_to_source"]-1)<1e-13
        assert abs(x["T_GeV"]/y["T_GeV"]-1)<1e-12
    return dict(scope="Matched saved states of the short independent warm continuation, "
                      "not a proof of pointwise occupations between these saved times, "
                      "nor of occupations in older omitted tails. Heavy fraction is a "
                      "frozen stationary diagnostic, not a coupled heavy trajectory.",
                cases=[dict(q=16,**a),dict(q=32,**b)])


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("q16_previous","q32_previous","q16_tape","q32_tape"):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    print(json.dumps(run(*vars(a).values()),indent=2,allow_nan=False))
