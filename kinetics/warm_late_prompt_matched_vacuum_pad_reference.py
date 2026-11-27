"""Evolve the late prompt channel with the same mesh but no exterior activation."""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from joint_prompt_two_moment_evolution_reference import Collision,Cells
from joint_resonant_centered_stiff_reference import advance
from joint_resonant_frw_split_reference import pack
from late_radiation_background_reference import MPL,TabulatedRadiation
from three_species_decay_onsager_reference import M,m,mp
from warm_compact_absolute_prompt_reference import LocalAbsolutePrompt
from warm_corrected_late_prompt_cut_aligned_reference import observable_weights,predictor
from warm_coupled_serial_restart_reference import atomic,sha
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom,pad


def velocity_moment(p,theta):
    out=[]
    for j,(mass,g) in enumerate(((mp,1),(m,4))):
        c=Cells(p,mass,g)
        k=c.low[:,None]+c.dk[:,None]*c.u
        E=k+mass
        y=c.f(theta[j],E,np.arange(len(p)-1)[:,None])
        w=c.w*c.vol[:,None]*y
        pp2=k*(k+2*mass)
        out.append([float(np.sum(w*pp2/E**2)/w.sum()),
                    float(np.sum(w*pp2)/w.sum()),float(np.sum(w*pp2**2)/w.sum())])
    return out


def run(eos,prior,source,matched,out,substeps):
    A=read(source)
    S=A["snapshot"]
    B=read(matched)
    G=B["snapshot"]
    assert A["snapshot_sha256"]==digest(S)
    assert B["hash"]==digest(G) and B["time"]==0.
    assert B["plan"]["artificial_empty_shell_filling"]==0.
    T=S["T_visible_GeV"]
    assert T==G["T_visible_GeV"] and S["scale_factor"]==G["scale_factor"]
    keys=("entropy_variable_scalar","entropy_variable_fermion")
    q=np.asarray(G["lab_momentum_edges_GeV"])
    p=np.asarray(S["lab_momentum_edges_GeV"])
    np.testing.assert_array_equal(p,q[:len(p)])
    orig=pad(S,q)
    U=pack(orig)
    rho=get_mom(q,[G[k] for k in keys])[2].sum()  # exactly the active test clock
    n0,_,E0=get_mom(q,orig)
    bg=TabulatedRadiation(eos,T)
    H=np.sqrt((bg.state(1.)[1]+rho)/3)/MPL
    cut=np.hypot(M,1./read(prior)["branches"][-1]["entropy_expansion_from_original_anchor"])
    op=LocalAbsolutePrompt(Collision(q),H,[],True)
    w=observable_weights(op,cut)
    plan=dict(eos_sha256=sha(eos),source_state_sha256=A["snapshot_sha256"],matching_active_grid_sha256=B["hash"],
              exact_matching_clock_from_active_case=True,substeps=substeps,driver_sha256=sha(Path(__file__)),
              compact_jacobian_sha256=sha(Path(__file__).with_name("warm_compact_absolute_prompt_reference.py")))
    out.mkdir(parents=True,exist_ok=True)
    track=[]
    times=np.r_[0.,np.geomspace(1e-11,1e-6,6)]
    for i,t in enumerate(times):
        dest=out/f"step-{i}.json"
        if dest.exists():
            R=read(dest)
            assert R["plan"]==plan and R["hash"]==digest(R["snapshot"]) and R["time"]==t
            U=pack([R["snapshot"][k] for k in keys])
        else:
            stat=[]
            if i:
                dt=(t-times[i-1])/substeps
                for _ in range(substeps):
                    try:
                        guess=predictor(op,U,dt)
                        U,s=advance(op,U,dt,T,initial_theta=guess,vacuum_scaled=True,balanced_merit=True)
                    except (RuntimeWarning,RuntimeError,AssertionError,ValueError,FloatingPointError) as problem:
                        U,s=advance(op,U,dt,T,reference_guess=True,vacuum_scaled=True,balanced_merit=True)
                        s["original_guess_failure"]=str(problem)
                    stat.append(s)
            shape=op.arrays(U)
            s=dict(S,lab_momentum_edges_GeV=q.tolist(),**{k:a.tolist() for k,a in zip(keys,shape,strict=True)})
            N,_,E=get_mom(q,shape)
            defect=[float(np.max(abs(N.sum(axis=1)/n0.sum(axis=1)-1.))),float(abs(E.sum()/E0.sum()-1.))]
            assert max(defect)<6e-12
            R=dict(plan=plan,hash=digest(s),snapshot=s,time=float(t),number_and_energy_defects=defect,
                   velocity2_and_p2_p4_by_species=velocity_moment(q,shape),
                   cut_excess_over_matched_clock_energy=float(w@op.moment(U)[0]/rho),
                   signed_cut_flux_over_matched_energy=float(w@op.collision(U)[0]/rho),
                   microsteps=stat)
            atomic(dest,R)
        track.append({k:v for k,v in R.items() if k!="snapshot"})
        print(json.dumps(dict(done=i,clock_grid_shells=op.n)),file=sys.stderr,flush=True)
    return dict(scope="Old cohorts remain gravitational only, with exactly the active run's fixed clock. "
                      "Same daughter mesh and no arbitrary vacuum floor. No other processes.",
                times=track)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","prior","source","matched","out"):
        p.add_argument(k,type=Path)
    p.add_argument("--substeps",type=int,default=3)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.prior,a.source,a.matched,a.out,a.substeps),indent=2,allow_nan=False))
