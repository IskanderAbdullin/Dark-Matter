"""Replay the actual turning interval with a locked clock and exact tail ledger."""

import argparse
import json
from pathlib import Path

import numpy as np

from joint_prompt_two_moment_evolution_reference import Collision
from joint_resonant_centered_stiff_reference import advance
from joint_resonant_frw_split_reference import pack
from late_radiation_background_reference import MPL,TabulatedRadiation
from three_species_decay_onsager_reference import M
from warm_compact_absolute_prompt_reference import LocalAbsolutePrompt
from warm_corrected_late_prompt_cut_aligned_reference import observable_weights,predictor
from warm_coupled_serial_restart_reference import atomic,sha
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_excess_energy_tail_reference import integrate
from warm_prompt_reentry_shape_audit_reference import get_mom


def run(eos,prior,seed,target,out,n):
    A=read(seed)
    S=A["snapshot"]
    assert A["hash"]==digest(S)
    assert A["plan"]["artificial_empty_shell_filling"]==0.
    lo=A["time"]
    assert target>lo and n>0
    H0=read(prior)["branches"][-1]
    cut=np.hypot(M,1./H0["entropy_expansion_from_original_anchor"])
    p=np.array(S["lab_momentum_edges_GeV"])
    T=S["T_visible_GeV"]
    keys=("entropy_variable_scalar","entropy_variable_fermion")
    z=[np.asarray(S[k]) for k in keys]
    N,K,E=get_mom(p,z)
    rho=E.sum()
    H=np.sqrt((TabulatedRadiation(eos,T).state(1.)[1]+rho)/3)/MPL
    op=LocalAbsolutePrompt(Collision(p),H,[],True)
    q=pack(z)
    weights=observable_weights(op,cut)
    plan=dict(eos_sha256=sha(eos),seed_sha256=A["hash"],target=target,locked_microsteps=n,
              source_sha256=sha(Path(__file__)),normalized_cut_GeV=cut,hard_GeV=A["plan"]["hard_GeV"])
    out.mkdir(parents=True,exist_ok=True)
    rows=[]
    h=(target-lo)/n
    for i,time in enumerate(np.linspace(lo,target,n+1)):
        path=out/f"step-{i}.json"
        if path.exists():
            R=read(path)
            assert R["plan"]==plan and R["time"]==time and R["hash"]==digest(R["snapshot"])
            q=pack([R["snapshot"][k] for k in keys])
        else:
            before=weights@op.moment(q)[0]/rho
            stat=None
            if i:
                reason=[]
                for method in ("moment predictor","own initial state","reference guess"):
                    try:
                        theta=(predictor(op,q,h) if method=="moment predictor" else
                               q.copy() if method=="own initial state" else None)
                        Z,stat=advance(op,q,h,T,initial_theta=theta,reference_guess=theta is None,
                                       vacuum_scaled=True,balanced_merit=True)
                        q=Z
                        break
                    except (RuntimeError,RuntimeWarning,ValueError,AssertionError,FloatingPointError) as e:
                        reason.append(str(e))
                else:
                    raise RuntimeError((i,reason))
                stat.update(seed_used=method,rejected_guesses=reason)
            D=dict(S,**{k:y.tolist() for k,y in zip(keys,op.arrays(q),strict=True)})
            W=weights@op.moment(q)[0]/rho
            f=weights@op.collision(q)[0]/rho
            delta=float(W-before-(h*f if i else 0.))
            assert abs(delta)<2e-10*max(abs(W),abs(before),1e-240)
            n1,k1,e1=get_mom(p,op.arrays(q))
            defects=[float(np.max(abs(n1.sum(axis=1)/N.sum(axis=1)-1.))),float(abs(e1.sum()/rho-1.))]
            assert max(defects)<2e-12
            independent=integrate(D,cut,8,True)/(rho*H)
            R=dict(plan=plan,hash=digest(D),snapshot=D,time=float(time),tail_energy=float(W),
                   signed_native_flux=float(f),signed_angular_flux=float(independent[3]),
                   angular_consumed_and_emitted_daughter_excess=independent[4:].tolist(),
                   exact_discrete_tail_balance_residual=delta,num_energy_defects=defects,accepted=stat)
            atomic(path,R)
        rows.append({k:v for k,v in R.items() if k!="snapshot"})
        print(json.dumps(dict(done=i,of=n,t=float(time))),flush=True,file=__import__("sys").stderr)
    pair=[]
    for a,b in zip(rows[:-1],rows[1:],strict=True):
        if a["signed_native_flux"]*b["signed_native_flux"]<0.:
            pair.append([a["time"],b["time"]])
    return dict(scope="Same cut-aligned phase and exact locked backward-Euler steps on a selected interval. "
                      "No cosmological heating, untracked exterior, or uniqueness of the first turning point.",
                sign_change_brackets=pair,rows=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","prior","seed","out"):
        p.add_argument(k,type=Path)
    p.add_argument("--target",type=float,default=1e-7)
    p.add_argument("--substeps",type=int,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.prior,a.seed,a.target,a.out,a.substeps),indent=2,allow_nan=False))
