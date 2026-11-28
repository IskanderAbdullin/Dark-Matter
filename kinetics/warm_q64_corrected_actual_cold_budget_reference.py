"""Replay a corrected guided q64 interval and close the actual cold-cell budget."""

import argparse
import json
from pathlib import Path

import numpy as np

import coupled_warm_pair_continuation_reference as kernel
from compact_prompt_jacobian_reference import local_prompt
from joint_resonant_centered_stiff_reference import advance
from joint_resonant_frw_split_reference import pack
from warm_cold_light_pair_partition_reference import split_rate
from warm_coupled_serial_restart_reference import sha
from warm_light_occupation_stage_audit_reference import digest, read


def run(eos,start,end,narrow):
    A,B,G=map(read,(start,end,narrow))
    for R in (A,B,G):
        assert digest(R["snapshot"])==R["snapshot_sha256"]
    assert B["parent_snapshot_sha256"]==A["snapshot_sha256"]
    assert B["plan"]==A["plan"] and B["plan"]["eos_sha256"]==sha(eos)
    assert B["plan"]["guide_snapshots"][B["step"]-1]==G["snapshot_sha256"]
    assert B["audit"]["trial"]["used_q32_guess"]
    S,U,K=[R["snapshot"] for R in (A,B,G)]
    target=[np.asarray(K["entropy_variable_"+k]) for k in ("scalar","fermion")]
    b=U["scale_factor"]/S["scale_factor"]
    dt=np.log(b)
    pair, drift = [], []
    stuff={}
    old_rate,old_collide,old_redshift=kernel.rate,kernel.collide,kernel.redshift

    def volume(p):
        return 4*p[1]**3/(6*np.pi**2)

    def occ(p,th):
        return kernel.moments(p,th)[1,0,0]/volume(p)

    def redshift(p,b0,b1,u):
        R=old_redshift(p,b0,b1,u)
        drift.append(float(occ(p/b1,R[0])-occ(p/b0,u)))
        return R

    def rate(p,theta,native):
        r,e,d,X=split_rate(p,theta,native,(4,8,6,8,8),cut_cells=9)
        pair.append(float(r[1,0,1]/volume(p)))
        return r,e,d

    def collide(model,u,h,T):
        assert model.prompt and model.local is not None
        initial=[w.copy() for w in model.arrays(u)]
        for a,z in zip(initial,target,strict=True):
            a[:,:z.shape[1]]=z
        detail={}
        Z,res=advance(model,u,h,T,reference_guess=False,vacuum_scaled=True,
                      initial_theta=pack(initial),balanced_merit=True,centered_diagnostic=detail)
        ref,q,r=(detail[k] for k in ("reference","centered","parent_ratio"))
        E,_=model.elastic(q)
        R,_=local_prompt(model,q,ref,r,model.local)
        mixed,_=model.centered(q,ref,r)
        w=2*model.n
        V=model.k.c[1].vol[0]
        np.testing.assert_allclose(mixed[w],E[w]+R[w],rtol=2e-14,atol=0.)
        before=model.moment(u)[0][w]/V
        after=model.moment(Z)[0][w]/V
        stuff.update(dt_hubble=h,H_GeV=model.H,volume=V,stiff_actual_change=after-before,
                     direct_prompt_change=h*R[w]/V,elastic_change=h*E[w]/V,
                     stiff_count_defect=after-before-h*(E[w]+R[w])/V,
                     residual=res["residual"])
        return Z,res["residual"],1,0

    after={}
    try:
        kernel.rate,kernel.collide,kernel.redshift=rate,collide,redshift
        report=kernel.run(eos,S,1,b,checkpoint=after,track_exterior=True,
                          local_prompt_jacobian=True,scheduled_single_block=True)
    finally:
        kernel.rate,kernel.collide,kernel.redshift=old_rate,old_collide,old_redshift
    assert len(pair)==4 and len(drift)==2
    p=np.asarray(after["lab_momentum_edges_GeV"])
    cur=[np.asarray(after["entropy_variable_"+k]) for k in ("scalar","fermion")]
    ref=[np.asarray(U["entropy_variable_"+k]) for k in ("scalar","fermion")]
    nx,ny=kernel.moments(p,cur),kernel.moments(p,ref)
    shell_l1=abs(nx-ny).sum(axis=-1)/ny.sum(axis=-1)
    assert np.max(shell_l1)<2e-11
    old=[np.asarray(S["entropy_variable_"+k]) for k in ("scalar","fermion")]
    exact=occ(p,cur)-occ(np.asarray(S["lab_momentum_edges_GeV"]),old)
    kicks=.5*dt*np.array(pair)[1::2]/stuff["H_GeV"]
    residual=exact-kicks.sum()-sum(drift)-stuff["stiff_actual_change"]
    assert abs(float(residual))<2e-16
    return dict(scope="Exact consecutive selected corrected q64 replay and discrete operator budget, "
                      "not a continuum sign or an omitted-channel bound.",
                start_sha256=A["snapshot_sha256"],end_sha256=B["snapshot_sha256"],
                narrow_guess_sha256=G["snapshot_sha256"],replay_sha256=digest(after),
                checkpoint_temperature_GeV=after["T_visible_GeV"],
                maximum_archived_endpoint_shell_L1=float(np.max(shell_l1)),
                accepted_energy_and_count_ledger=[report["energy_ledger"],report["count_ledger"]],
                redshift_changes=drift,pair_kick_changes=kicks.tolist(),stiff=stuff,
                observed_total_change=float(exact),unattributed_roundoff=float(residual),
                per_lna=dict(total=float(exact/dt),pair=float(sum(kicks)/dt),
                             prompt=float(stuff["direct_prompt_change"]/dt),
                             below_pole_elastic=float(stuff["elastic_change"]/dt)))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for x in ("eos","start","endpoint","q32_guide"):
        p.add_argument(x,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.start,a.endpoint,a.q32_guide),indent=2,allow_nan=False))
