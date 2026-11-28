"""Follow two archive-moment-matched shapes near the frozen current reversal."""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from joint_prompt_two_moment_evolution_reference import Collision
from joint_resonant_centered_stiff_reference import advance
from joint_resonant_frw_split_reference import pack
from late_radiation_background_reference import MPL,TabulatedRadiation
from three_species_decay_onsager_reference import M
from warm_compact_absolute_prompt_reference import LocalAbsolutePrompt
from warm_corrected_exterior_reactivated_tail_stress_reference import reactivate
from warm_corrected_late_prompt_cut_aligned_reference import excess,observable_weights,predictor,split_at_cut
from warm_coupled_serial_restart_reference import atomic,sha
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_excess_energy_tail_reference import integrate
from warm_prompt_reentry_shape_audit_reference import get_mom


def run(eos,prior,archive,out,hard,n,halfwidth):
    Y=read(archive)
    S=Y["snapshot"]
    assert Y["snapshot_sha256"]==digest(S)
    C=np.hypot(M,1./read(prior)["branches"][-1]["entropy_expansion_from_original_anchor"])
    A,_,_,_,_=reactivate(S,hard,halfwidth,cap=1.05)
    A=split_at_cut(A,C)
    p=np.asarray(A["lab_momentum_edges_GeV"])
    T=A["T_visible_GeV"]
    keys=("entropy_variable_scalar","entropy_variable_fermion")
    th=[np.asarray(A[k]) for k in keys]
    N0,_,E0=get_mom(p,th)
    rho=E0.sum()
    H=np.sqrt((TabulatedRadiation(eos,T).state(1.)[1]+rho)/3.)/MPL
    u=pack(th)
    op=LocalAbsolutePrompt(Collision(p),H,[],True,stream_local=True)
    weights=observable_weights(op,C)
    start_W=float(weights@op.moment(u)[0]/rho)
    lineage=dict(archive_sha256=Y["snapshot_sha256"],seed_sha256=digest(A),eos_sha256=sha(eos),
                 normalized_cut_GeV=float(C),hard_GeV=hard,shell_halfwidth_GeV=halfwidth,
                 substeps=n,solver_source_sha256=sha(Path(__file__)))
    rows=[]
    clock=[0.,1e-11,1e-10,1e-9,1e-8,3e-8,1e-7]
    out.mkdir(parents=True,exist_ok=True)
    def write(t,diag):
        P=dict(A,**{k:x.tolist() for k,x in zip(keys,op.arrays(u),strict=True)})
        N,_,E=get_mom(p,op.arrays(u))
        numdef=float(max(abs(N.sum(axis=1)/N0.sum(axis=1)-1.)))
        edef=float(abs(E.sum()/rho-1.))
        assert max(numdef,edef)<8e-12
        w=float(weights@op.moment(u)[0]/rho)
        w_direct=sum(excess(P,C))/rho
        assert abs(w-w_direct)<4e-11*max(w,w_direct,1e-70)
        prompt=float(weights@op.collision(u)[0]/rho)
        angle=(integrate(P,C,8,True)/(rho*H)).tolist()
        assert abs(angle[3]-prompt)<3e-5*max(angle[1:3])
        row=dict(time=t,plan=lineage,hash=digest(P),tail_energy=w,signed_native_flux=prompt,
                 angular_signed_flux=angle[3],discrete_N_and_E_errors=[numdef,edef],substep_solver=diag)
        atomic(out/f"step-{len(rows)}.json",dict(row,snapshot=P))
        rows.append(row)
    write(0.,[])
    for a,b in zip(clock,clock[1:]):
        files=[]
        for _ in range(n):
            h=(b-a)/n
            old_w=weights@op.moment(u)[0]/rho
            guide=predictor(op,u,h)
            u,info=advance(op,u,h,T,initial_theta=guide,vacuum_scaled=True,balanced_merit=True)
            new=weights@op.moment(u)[0]/rho
            flux=weights@op.collision(u)[0]/rho
            defect=abs(new-old_w-h*flux)
            assert defect<4e-10*max(new,old_w,1e-70)
            info["projected_W_balance_error"]=float(defect)
            files.append(info)
        write(b,files)
        print(hard,n,b,rows[-1]["signed_native_flux"],file=sys.stderr,flush=True)
    atomic(out/"summary.json",dict(plan=lineage,W_initial=start_W,H_GeV=float(H),rows=rows,
                                  scope="Fixed late phase; prompt only, stationary parent, closed 1.05 GeV support."))


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","prior","archive","out"):
        ap.add_argument(k,type=Path)
    ap.add_argument("--hard",type=float,required=True)
    ap.add_argument("--substeps",type=int,default=3)
    ap.add_argument("--halfwidth",type=float,default=.005)
    a=ap.parse_args()
    assert a.substeps>0 and a.halfwidth>0.
    run(a.eos,a.prior,a.archive,a.out,a.hard,a.substeps,a.halfwidth)
