"""Run or audit a matched coupled q32/q64 support control."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import expit

from frw_exterior_restart_reference import import_restart
from joint_prompt_two_moment_evolution_reference import Cells
from three_species_decay_onsager_reference import m,mp
from warm_coupled_serial_restart_reference import bytes_of,main,sha
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom,pad


def sources(folder,D):
    a=read(D/"warm-light-occupation-next-q32-tape-2026-11-26.json.gz")[0]
    A=a["snapshot"]
    assert digest(A)==a["snapshot_sha256"]
    p=np.asarray(A["lab_momentum_edges_GeV"])
    q=np.r_[p,p[-1]*np.geomspace(1.,2.,6)[1:]]
    t=pad(A,q)
    B=dict(A,lab_momentum_edges_GeV=q.tolist(),entropy_variable_scalar=t[0].tolist(),
           entropy_variable_fermion=t[1].tolist())
    folder.mkdir(parents=True,exist_ok=True)
    for name,S in (("source32.json",A),("source64.json",B)):
        path=folder/name
        b=bytes_of({"snapshots":[S]})
        if path.exists() and path.read_bytes()!=b:
            raise ValueError("Mismatched source or destination.")
        if not path.exists():
            path.write_bytes(b)
    return A,B,p,q


def audit(folder,D,eos,nsteps):
    A,B,p,q=sources(folder,D)
    ix=len(p)-1
    a0=get_mom(p,pad(A,p))[0].sum(axis=1)*A["scale_factor"]**3
    x0=get_mom(p,pad(A,p))[2].sum()*A["scale_factor"]**3
    C=import_restart(A["positive_exterior_restart"]).array()
    results=[]; arrays=[]; spectra=[]
    for n,src in ((32,A),(64,B)):
        f=folder/f"q{n}"/"checkpoint.json"
        k=read(f)
        assert k["blocks_finished"]==k["blocks_total"]==nsteps and len(k["ledger"])==nsteps
        assert k["plan"]["source_sha256"]==sha(folder/f"source{n}.json")
        assert k["plan"]["eos_sha256"]==sha(eos)
        assert k["plan_sha256"]==digest(k["plan"])
        s=k["snapshot"]
        assert digest(s)==k["snapshot_sha256"]
        N,K,E=get_mom(np.asarray(s["lab_momentum_edges_GeV"]),
                      [s["entropy_variable_scalar"],s["entropy_variable_fermion"]])
        a=s["scale_factor"]
        b=a/src["scale_factor"]
        assert abs(b/1.001-1)<5e-13
        ext=import_restart(s["positive_exterior_restart"]).array()
        core=len(C.T)
        save=ext[:,:core].copy()
        save[[0,5]]*=np.array([b,b**3])[:,None]
        np.testing.assert_allclose(save,C,rtol=3e-14,atol=0.)
        extra=ext[:,core:]
        rmax=max(max(abs(u["energy_ledger"]),abs(u["count_ledger"])) for u in k["ledger"])
        assert rmax<1e-10
        inclusive=np.array([np.sum(extra[5][extra[1]==x])*a**3 for x in (mp,m)])
        total=(N.sum(axis=1)*a**3+inclusive).sum()
        assert abs(total/a0.sum()-1)<1e-10
        mass=np.array([mp,m])
        new_exterior=(extra[5]*extra[6])*a**3  # birth energies, not endpoint energy
        pressure=[]; spectrum=[]
        for j,d in enumerate(mass):
            mesh=np.asarray(s["lab_momentum_edges_GeV"])
            c=Cells(mesh,d,(1,4)[j])
            t=np.asarray(s[("entropy_variable_scalar","entropy_variable_fermion")[j]])
            x=t[0,:,None]+t[1,:,None]*c.u
            ff=expit(-x) if j else np.exp(-x)/(-np.expm1(-x))
            z=t[0,:ix,None]+t[1,:ix,None]*np.linspace(0.,1.,17)
            spectrum.append(expit(-z) if j else np.exp(-z)/(-np.expm1(-z)))
            e=d+c.low[:,None]+c.dk[:,None]*c.u
            Q=(e-d)*(e+d)
            pressure.append(((c.w*ff*(Q/(3*e))).sum(axis=1)*c.vol*a**3)[:ix])
        spectra.append(np.asarray(spectrum))
        arrays.append((N[:,:ix]*a**3,E[:,:ix]*a**3,np.asarray(pressure)))
        results.append(dict(q=n,checkpoint_sha256=k["snapshot_sha256"],T_GeV=s["T_visible_GeV"],
                            scale_factor=a,latest_block_ledger_defect=rmax,
                            inclusive_combined_number_relative_defect=abs(total/a0.sum()-1),
                            extra_represented_beyond_q32_number_relative_initial=(N[:,ix:].sum(axis=1)*a**3/a0).tolist(),
                            extra_represented_beyond_q32_energy_relative_initial=(E[:,ix:].sum(axis=1)*a**3/x0).tolist(),
                            new_collisionless_receipt_count=extra.shape[1],
                            new_collisionless_receipt_species_number_relative_initial=(inclusive/a0).tolist(),
                            new_collisionless_birth_energy_relative_initial=float(new_exterior.sum()/x0)))
    assert abs(results[0]["T_GeV"]/results[1]["T_GeV"]-1)<3e-14
    rel=spectra[1]-spectra[0]
    return dict(shared_max_occupation_differences_sampled=np.max(abs(rel),axis=(1,2)).tolist(),
                first_light_shell_occupation_at_rest=[float(s[1,0,0]) for s in spectra],
                scope="Exact inherited collisionless receipt gravity shared in both branches. Vacuum "
                      "explicit q64 padding is not a reconstruction of older exterior population; "
                      "only new represented populations react. Coupled prompt, selected elastic, "
                      "and signed scalar pair, not complete cosmological evolution.",
                blocks=nsteps,controls=results,shared_q32_number_L1_over_initial=(np.sum(abs(arrays[1][0]-arrays[0][0]),axis=1)/a0).tolist(),
                shared_q32_energy_L1_over_initial=(np.sum(abs(arrays[1][1]-arrays[0][1]),axis=1)/x0).tolist(),
                shared_q32_pressure_L1_over_initial_energy=(np.sum(abs(arrays[1][2]-arrays[0][2]),axis=1)/x0).tolist(),
                within_shared_binwise=dict(number=[(s[0]/a0[:,None]).tolist() for s in arrays],
                                           energy=[(s[1]/x0).tolist() for s in arrays],
                                           pressure=[(s[2]/x0).tolist() for s in arrays]))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("work",type=Path)
    p.add_argument("--tapes",type=Path,default=Path("results"))
    p.add_argument("--steps",type=int,choices=(2,4,8),default=2)
    p.add_argument("--refined-work",type=Path)
    p.add_argument("--execute",action="store_true")
    a=p.parse_args()
    if a.execute:
        sources(a.work,a.tapes)
        for q in (32,64):
            main(a.eos,a.work/f"source{q}.json",a.work/f"q{q}",1.001,a.steps,True,None)
    base=audit(a.work,a.tapes,a.eos,a.steps)
    if a.refined_work is not None:
        assert a.steps==2
        fine=audit(a.refined_work,a.tapes,a.eos,4)
        time={}
        for key,val in fine["within_shared_binwise"].items():
            x=np.asarray(base["within_shared_binwise"][key])
            y=np.asarray(val)
            time[key]=np.sum(abs(x-y),axis=2).tolist()
        base=dict(two_block=base,four_block=fine,
                  time_discretization_L1_by_q32_q64=dict(**time))
    print(json.dumps(base,indent=2,allow_nan=False))
