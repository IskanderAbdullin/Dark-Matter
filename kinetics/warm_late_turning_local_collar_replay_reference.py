"""Refine only the threshold collar and replay the exact selected clock."""

import argparse
import json
from pathlib import Path

import numpy as np

from joint_prompt_two_moment_evolution_reference import Collision,kin
from joint_resonant_centered_stiff_reference import advance
from joint_resonant_frw_split_reference import pack
from late_radiation_background_reference import MPL,TabulatedRadiation
from three_species_decay_onsager_reference import M,m,mp
from warm_compact_absolute_prompt_reference import LocalAbsolutePrompt
from warm_corrected_late_prompt_cut_aligned_reference import observable_weights,predictor
from warm_coupled_serial_restart_reference import atomic,sha
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom


def remap(p,q,S):
    i=np.searchsorted(p,(q[:-1]+q[1:])/2,side="right")-1
    assert np.all(p[i]<=q[:-1]) and np.all(q[1:]<=p[i+1])
    out=[]
    for d,k in ((mp,"scalar"),(m,"fermion")):
        x=kin(p,d);y=kin(q,d);a,b=np.asarray(S["entropy_variable_"+k])
        w=np.diff(x)[i]
        out.append(np.array([a[i]+b[i]*(y[:-1]-x[i])/w,b[i]*np.diff(y)/w]))
    return out


def run(eos,trail,out):
    A=read(trail/"step-0.json")
    S=A["snapshot"]; assert A["hash"]==digest(S)
    p=np.asarray(S["lab_momentum_edges_GeV"])
    select=(p[:-1]>=.645)&(p[1:]<=.765)
    assert 2<=select.sum()<=9
    q=np.sort(np.r_[p,(p[:-1]+p[1:])[select]/2])
    base=remap(p,q,S)
    U=pack(base)
    N,K,E=get_mom(q,base)
    rho=E.sum()
    H=np.sqrt((TabulatedRadiation(eos,S["T_visible_GeV"]).state(1.)[1]+rho)/3)/MPL
    op=LocalAbsolutePrompt(Collision(q),H,[],True)
    cut=A["plan"]["normalized_cut_GeV"]; weights=observable_weights(op,cut)
    plan=dict(eos_sha256=sha(eos),input_sha256=A["hash"],cut_GeV=cut,driver_sha256=sha(Path(__file__)))
    out.mkdir(parents=True,exist_ok=True)
    start=A; rows=[]
    for j in range(1,7):
        end=read(trail/f"step-{j}.json")
        assert end["hash"]==digest(end["snapshot"]) and end["plan"]==A["plan"]
        h=end["time"]-start["time"]
        seed=pack(remap(p,q,end["snapshot"]))
        old=weights@op.moment(U)[0]/rho
        try:
            v,stat=advance(op,U,h,S["T_visible_GeV"],initial_theta=seed,vacuum_scaled=True,balanced_merit=True)
        except (ValueError,AssertionError,RuntimeError,RuntimeWarning,FloatingPointError):
            v,stat=advance(op,U,h,S["T_visible_GeV"],initial_theta=predictor(op,U,h),
                           vacuum_scaled=True,balanced_merit=True)
        U=v
        th=op.arrays(U)
        phase=dict(S,lab_momentum_edges_GeV=q.tolist(),entropy_variable_scalar=th[0].tolist(),entropy_variable_fermion=th[1].tolist())
        now=weights@op.moment(U)[0]/rho
        rate=weights@op.collision(U)[0]/rho
        defect=abs(now-old-h*rate)
        assert defect<2e-10*max(old,now,1e-200)
        row=dict(plan=plan,hash=digest(phase),snapshot=phase,time=end["time"],tail_excess_over_rho=float(now),
                 native_signed_flux=float(rate),guide_was_initial_guess_sha256=end["hash"],solver=stat,
                 same_clock_original_tail_excess=end["tail_energy"],discrete_tail_balance_residual=float(defect))
        atomic(out/f"step-{j}.json",row)
        print(json.dumps(dict(done=j,time=end["time"])),file=__import__("sys").stderr,flush=True)
        rows.append({k:v for k,v in row.items() if k!="snapshot"})
        start=end
    return dict(scope="Targeted collar split only. Original endpoints are Newton guesses, not evolved inputs.",
                old_and_split_shell_counts=[len(p)-1,len(q)-1],rows=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","trail","output"):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.trail,a.output),indent=2,allow_nan=False))
