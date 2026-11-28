"""Conditional reactivation of known cohort moments on the matched q64 mesh."""

import argparse
import json
from pathlib import Path

import numpy as np

from frw_exterior_cohort_reference import Exterior
from frw_exterior_restart_reference import export_restart,import_restart
from joint_prompt_two_moment_evolution_reference import Cells
from three_species_decay_onsager_reference import m,mp
from warm_coupled_serial_restart_reference import bytes_of,main,sha
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom,pad


def seed(s):
    p=np.asarray(s["lab_momentum_edges_GeV"])
    ix=len(p)-1
    wide=np.r_[p,p[-1]*np.geomspace(1.,2.,6)[1:]]
    th=pad(s,wide)
    ext=import_restart(s["positive_exterior_restart"])
    Q=ext.array()
    a,mm,_,_,_,n,e=Q
    now=a*np.sqrt((e-mm)*(e+mm))
    energy=np.hypot(now,mm)
    selected=(now>=p[-1])&(now<wide[-1])
    klass=np.searchsorted(wide,now,side="right")-1
    assert np.any(selected)
    deltas=[]
    for j,d in enumerate((mp,m)):
        ok=selected&(mm==d)
        N=np.bincount(klass[ok],weights=n[ok],minlength=len(wide)-1)
        K=np.bincount(klass[ok],weights=n[ok]*now[ok]**2/(energy[ok]+d),minlength=len(wide)-1)
        c=Cells(wide,d,(1,4)[j])
        out,err,_=c.reconstruct([N,K])
        assert err<5e-10
        assert not np.any(N[:ix])
        th[j][:,ix:]=np.asarray(out)[:,ix:]
        deltas.append([N.sum(),K.sum()+d*N.sum()])
    inactive=Exterior()
    inactive.cohorts=[x for x,kept in zip(ext.cohorts,selected,strict=True) if not kept]
    inactive.injected=float((a[~selected]*n[~selected])@e[~selected])
    inactive.trace=max(0.,float(n[~selected]@energy[~selected]-inactive.injected))
    assert abs(inactive.energy(1.)+np.sum(np.asarray(deltas)[:,1])-ext.energy(1.))<1e-12*ext.energy(1.)
    N,_,E=get_mom(wide,th)
    T=E.sum()+inactive.energy(1.)
    old_mom=get_mom(p,pad(s,p))
    P=old_mom[2].sum()+ext.energy(1.)
    before=old_mom[0].sum(axis=1)+np.array([n[mm==d].sum() for d in (mp,m)])
    after=N.sum(axis=1)+np.array([n[(mm==d)&(~selected)].sum() for d in (mp,m)])
    assert np.max(abs(after/before-1))<2e-12
    assert abs(T/P-1)<2e-12
    S=dict(s,lab_momentum_edges_GeV=wide.tolist(),entropy_variable_scalar=th[0].tolist(),
           entropy_variable_fermion=th[1].tolist(),
           represented_species_count_a3=(N.sum(axis=1)*s["scale_factor"]**3).tolist(),
           positive_exterior_restart=export_restart(inactive))
    return S,dict(activated_cohorts=int(selected.sum()),still_collisionless=int((~selected).sum()),
                  moment_closure_fractions=(np.asarray(deltas)/np.c_[N.sum(axis=1),np.repeat(T,2)]).tolist(),
                  total_energy_relative_transfer_defect=abs(T/P-1.),
                  total_species_number_transfer_defect=np.max(abs(after/before-1)).item())


def run(D,eos,work,baseline,execute):
    a=read(D/"warm-light-occupation-next-q32-tape-2026-11-26.json.gz")[0]
    s=a["snapshot"]
    assert a["snapshot_sha256"]==digest(s)
    S,ledger=seed(s)
    work.mkdir(parents=True,exist_ok=True)
    blob=bytes_of({"snapshots":[S]})
    data=work/"seed.json"
    if data.exists():
        assert data.read_bytes()==blob
    else:
        data.write_bytes(blob)
    if execute:
        main(eos,data,work/"q64",1.001,2,True,None)
    z=read(work/"q64"/"checkpoint.json")
    x=read(baseline)
    assert z["blocks_total"]==x["blocks_total"]==2 and z["blocks_finished"]==x["blocks_finished"]==2
    assert z["plan"]["source_sha256"]==sha(data)
    A,B=[v["snapshot"] for v in (z,x)]
    for q,checkpoint in ((A,z),(B,x)):
        assert digest(q)==checkpoint["snapshot_sha256"]
    assert abs(A["T_visible_GeV"]/B["T_visible_GeV"]-1)<1e-13
    idx=len(s["lab_momentum_edges_GeV"])-1
    N0,_,E0=get_mom(np.asarray(s["lab_momentum_edges_GeV"]),pad(s,np.asarray(s["lab_momentum_edges_GeV"])))
    M=[]
    for V in (A,B):
        p=np.asarray(V["lab_momentum_edges_GeV"])
        M.append(get_mom(p,[V["entropy_variable_scalar"],V["entropy_variable_fermion"]]))
    a3=(A["scale_factor"]/s["scale_factor"])**3
    scaler=N0.sum(axis=1)
    core=np.sum(abs(M[0][0][:,:idx]-M[1][0][:,:idx]),axis=1)*a3/scaler
    core_e=np.sum(abs(M[0][2][:,:idx]-M[1][2][:,:idx]),axis=1)*a3/E0.sum()
    totals=[(v[0][:,idx:].sum(axis=1)*a3/scaler).tolist() for v in M]
    return dict(scope="One representative momentum per earlier stored cohort, matched by its birth "
                      "number and mean energy and redshifted to the anchor. This realizes the "
                      "one-momentum Jensen lower energy, not the physical unresolved spread. "
                      "Selected cohorts are removed from the collisionless receipt ledger before "
                      "being permitted to react. Coupled selected-operator control only.",
                seed=ledger,new_seed_sha256=sha(data),reentry_checkpoint_sha256=z["snapshot_sha256"],
                vacuum_checkpoint_sha256=x["snapshot_sha256"],
                shared_core_counts_L1_over_initial=core.tolist(),
                shared_core_energy_L1_over_initial_total=core_e.tolist(),
                actual_tail_count_over_original=[dict(reactivated=totals[0],vacuum=totals[1])],
                largest_native_defect=max(max(abs(d["energy_ledger"]),abs(d["count_ledger"])) for d in z["ledger"]))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for lab in ("eos","work","baseline"):p.add_argument(lab,type=Path)
    p.add_argument("--tapes",type=Path,default=Path("results"))
    p.add_argument("--execute",action="store_true")
    x=p.parse_args()
    print(json.dumps(run(x.tapes,x.eos,x.work,x.baseline,x.execute),indent=2,allow_nan=False))
