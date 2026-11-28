"""Prompt and elastic dark evolution plus signed pair conversion after injection."""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

from frw_exterior_cohort_reference import Exterior
from frw_exterior_restart_reference import export_restart, import_restart
from joint_below_pole_elastic_linear_modes_reference import setup
from joint_prompt_two_moment_evolution_reference import Collision
from joint_resonant_frw_split_reference import pack,redshift
from late_archived_shape_continuation_reference import Late,collide,report
from late_radiation_background_reference import MPL,TabulatedRadiation
from scalar_pair_expanding_continuation_reference import MASS,e_total,moments,rate,reconstruct
from joint_nonlinear_prompt_return_moments_reference import P


def run(table,s,n,bmax,pair=True,checkpoint=None,progress_every=0,track_exterior=False,local_prompt_jacobian=False,scheduled_single_block=False):
    q=np.array(s["lab_momentum_edges_GeV"])
    t0=s["T_visible_GeV"]
    ar=s["scale_factor"]
    clock=TabulatedRadiation(table,t0)
    th=[np.asarray(s[k]) for k in ("entropy_variable_scalar","entropy_variable_fermion")]
    ei=moments(q,th)
    n0=ei[:,0].sum(axis=1)
    e0=e_total(ei)
    assert max(abs(n0*ar**3/np.asarray(s["represented_species_count_a3"])-1.))<3e-12
    prior=s.get("positive_exterior_restart") if track_exterior else None
    ex=import_restart(prior) if prior is not None else Exterior()
    e_ex0, trace_ex0=ex.energy(1.),ex.trace
    if track_exterior and not np.all(np.isin(ex.array()[1],MASS)):
        raise ValueError("Unexpected mass in predecessor exterior archive.")
    exn=np.zeros(2); pair_count=np.zeros(2)
    work=0.; peak=0.; mx=0.; errors=[];accepted=0;rejected=0
    grid=np.unique(np.r_[np.linspace(0.,np.log(bmax),n+1),np.clip(np.log(q[-1]/P),0.,np.log(bmax))])
    if scheduled_single_block:
        # The controlling outer grid already split any interior prompt crossing.
        # At an endpoint only, suppress a duplicate from a few ulps of rebasing.
        assert n==1
        cut=np.log(q[-1]/P)
        assert cut<=2e-14 or cut>=np.log(bmax)-2e-14
        grid=np.array([0.,np.log(bmax)])
    rows=[]
    for x,y in zip(grid[:-1],grid[1:]):
        a,b,c=np.exp([x,(x+y)/2,y])
        t=clock.state(b)[0]
        th,wr,*_=redshift(q,a,b,th);work+=wr.sum()
        p=q/b
        cur=moments(p,th)
        ex.advance(b)
        H=np.sqrt((clock.state(b)[1]+e_total(cur)+ex.energy(b)/b**4)/3)/MPL
        half=(y-x)/H/2
        def kick(state):
            nonlocal peak
            u=moments(p,state)
            rr,_,_=rate(p,state,(3,7,6,2))
            trial=reconstruct(p,u+.5*half*rr[:,:,1:-1])
            rr,err,z=rate(p,trial,(3,7,6,2))
            errors.append(err)
            peak=max(peak,z)
            delta=half*rr
            d=b**3*delta[:,0,1:-1].sum(axis=1)
            pair_count[:]+=d
            exn[:]+=b**3*delta[:,0,[0,-1]].sum(axis=1)
            for j,mi in enumerate(MASS):
                for side,k in enumerate((0,-1)):
                    ex.add(b,mi,3+j if track_exterior else j,side,p[k],b**3*delta[j,0,k],b**4*delta[j,1,k])
            return reconstruct(p,u+delta[:,:,1:-1])
        if pair:
            th=kick(th)
        cell=Collision(p)
        prompt=b<q[-1]/P and len(cell.y)>0
        seed=Late(cell,H,[],prompt)
        u=pack(th)
        _,_,T=seed.reference(u,t)
        panels=[]
        setup(max(T,t),cell,T,seed.eta,3,2,2,8,32.,.82,cache=panels,cache_only=True)
        model=Late(cell,H,panels,prompt,compact=True,local_jac=local_prompt_jacobian)
        next_state,res,ok,no=collide(model,u,y-x,T)
        mx=max(mx,res);accepted+=ok;rejected+=no
        th=model.arrays(next_state)
        if pair:
            th=kick(th)
        th,wr,*_=redshift(q,b,c,th);work+=wr.sum()
        ex.advance(c)
        out=moments(q/c,th)
        nf=out[:,0].sum(axis=1)*c**3
        r,_=report(Late(Collision(q/c),H,[],False),th,c,n0+pair_count)
        r["retained_counts_relative_to_start"]=(nf/n0).tolist()
        r["inclusive_counts_relative_to_start"]=((nf+exn)/n0).tolist()
        r["cold_count_relative_to_initial_total_scalar"]=out[0,0,0]*c**3/n0[0]
        r["photon_temperature_GeV"]=clock.state(c)[0]
        rows.append(r)
        j=len(rows)
        if progress_every and (j%progress_every==0 or j==len(grid)-1):
            progress=dict(kind="warm_pair_checkpoint",steps_requested=n,
                          blocks_finished=j,blocks_total=len(grid)-1,
                          relative_expansion=float(c),photon_temperature_GeV=r["photon_temperature_GeV"],
                          stiff_microsteps_accepted=accepted,stiff_microsteps_rejected=rejected)
            print(json.dumps(progress,allow_nan=False),file=sys.stderr,flush=True)
    exr=ex.result()
    rel_e=(e_total(out)*c**4+(exr["exterior_a4_energy_one_mean_per_cohort"]-e_ex0)
           -work-(ex.trace-trace_ex0)-e0)/e0
    rel_n=(nf+exn).sum()/n0.sum()-1.
    exchange=(pair_count+exn).sum()/n0.sum()
    assert max(abs(rel_e),abs(rel_n),abs(exchange),max(r["count_defect"] for r in rows))<5e-9
    if checkpoint is not None:
        checkpoint.update(T_initial_GeV=s["T_initial_GeV"],T_visible_GeV=clock.state(c)[0],scale_factor=ar*c,
                lab_momentum_edges_GeV=(q/c).tolist(),entropy_variable_scalar=np.asarray(th[0]).tolist(),
                entropy_variable_fermion=np.asarray(th[1]).tolist(),represented_species_count_a3=(nf*ar**3).tolist(),
                steps=n)
        if track_exterior:
            checkpoint["positive_exterior_restart"]=export_restart(ex)
    if track_exterior:
        exr.update(predecessor_cohort_count=len(prior["cohorts"]) if prior else 0,
                   initial_relative_a4_exterior_energy=e_ex0,
                   increment_in_a4_exterior_trace=ex.trace-trace_ex0)
    scope="Radiation plus retained and newly produced exterior gravity; prompt and below-pole elastic medium plus optionally pair conversion. Prior exterior cohorts, visible slow source, heavy recycling and number-changing scalar self-interactions are not evolved."
    if track_exterior:
        scope=("The previous positive collisionless exterior contributes to expansion when available. "
               "It never enters tail collision occupations, and entry counts are not counted again. "
               "Visible slow source, heavy recycling and number-changing scalar self-interactions are not evolved.")
    return dict(steps=n,initial_temperature_GeV=t0,with_pair=pair,acceptances=accepted,rejections=rejected,
            maximum_stiff_residual=mx,energy_ledger=rel_e,count_ledger=rel_n,collision_number_ledger=exchange,
            max_atomic_pair_error=max(errors,default=0.),maximum_midpoint_inverse_over_forward=peak,
            exterior=exr,endpoint=rows[-1],trajectory=rows,
            scope=scope)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("table");p.add_argument("shape")
    p.add_argument("--stop",type=float,default=1.2)
    p.add_argument("--steps",type=int,nargs="+",default=[4,8])
    p.add_argument("--no-pair",action="store_true")
    p.add_argument("--save-snapshots")
    p.add_argument("--progress-every",type=int,default=0,
                   help="Write status after this many accepted outer blocks to stderr, without changing outputs.")
    p.add_argument("--local-prompt-jacobian",action="store_true",
                   help="Assemble the same prompt derivative in parent-local blocks; useful as a numerical control.")
    p.add_argument("--track-exterior",action="store_true",
                   help="Export/reuse positive collisionless exterior cohort gravity; not tail interactions.")
    a=p.parse_args()
    if a.track_exterior and not a.save_snapshots:
        p.error("--track-exterior requires --save-snapshots")
    assert a.stop>1 and min(a.steps)>0 and a.progress_every>=0
    s=json.loads(Path(a.shape).read_text())["snapshots"][-1]
    ck=[{} for _ in a.steps] if a.save_snapshots else [None for _ in a.steps]
    rows=[run(a.table,s,n,a.stop,not a.no_pair,checkpoint=k,progress_every=a.progress_every,track_exterior=a.track_exterior,
              local_prompt_jacobian=a.local_prompt_jacobian)
          for n,k in zip(a.steps,ck)]
    if a.save_snapshots:
        Path(a.save_snapshots).write_text(json.dumps(dict(snapshots=ck),indent=2,allow_nan=False)+"\n")
    print(json.dumps(dict(rows=rows),indent=2,allow_nan=False))
