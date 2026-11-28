"""Compare stiff integrators on an identical frozen prompt-plus-elastic shape."""
import argparse,json
from pathlib import Path
import hashlib

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

import numpy as np
from joint_below_pole_elastic_linear_modes_reference import setup
from joint_prompt_two_moment_evolution_reference import Collision
from joint_resonant_centered_stiff_reference import advance
from joint_resonant_frw_split_reference import pack
from late_archived_shape_continuation_reference import Late
from late_radiation_background_reference import MPL,TabulatedRadiation
from scalar_pair_expanding_continuation_reference import e_total,moments,reconstruct


def implicit(model,q,h,T):
    errors=[]
    for reference in (False,True):
        try:
            return advance(model,q,h,T,vacuum_scaled=True,reference_guess=reference)
        except (RuntimeError,FloatingPointError,RuntimeWarning,ValueError) as e:
            errors.append(e)
    raise errors[-1]


def midpoint(model,q,step,T,p):
    mid,diag=implicit(model,q,.5*step,T)
    U=moments(p,model.arrays(q))
    Mid=moments(p,model.arrays(mid))
    goal=2*Mid-U
    if np.any(goal[:,0]<=0.) or np.any(goal[:,1]<=0.):
        raise ValueError("Implicit midpoint reflection leaves the positive-moment cone.")
    out=pack(reconstruct(p,goal))
    actual=moments(p,model.arrays(out))
    scale=np.maximum(abs(actual),abs(goal))
    err=float(np.max(abs(goal-actual)/np.maximum(scale,1e-290)))
    _,S0,_=model.moment(q)
    _,S1,_=model.moment(out)
    if S1<S0*(1.-3e-13):
        raise ValueError("Implicit midpoint reflection decreased total entropy.")
    assert err<2e-9
    return out,diag


def sdirk2(model,q,h,T):
    g=1.-1./np.sqrt(2.)
    context={}
    U0=model.moment(q)[0]
    q1,d1=advance(model,q,g*h,T,vacuum_scaled=True,centered_diagnostic=context)
    F1,_=model.centered(context["centered"],context["reference"],context["parent_ratio"])
    q2,d2=advance(model,q,g*h,T,initial_theta=q1,vacuum_scaled=True,
                         target_moments=U0+(1.-g)*h*F1)
    d2["newton_iterations"]+=d1["newton_iterations"]
    d2["residual"]=max(d1["residual"],d2["residual"])
    return q2,d2

def run(eos,shape,step,N):
    assert tuple(N)==(2,4,8) and step>0
    s=json.loads(shape.read_text())["snapshots"][-1]
    p=np.asarray(s["lab_momentum_edges_GeV"])
    th=[np.asarray(s[k]) for k in ("entropy_variable_scalar","entropy_variable_fermion")]
    T0=s["T_visible_GeV"]
    H=np.sqrt((TabulatedRadiation(eos,T0).state(1.)[1]+e_total(moments(p,th)))/3.)/MPL
    k=Collision(p)
    assert len(k.y)>0
    u=pack(th)
    panels=[]
    seed=Late(k,H,[],True)
    _,_,T=seed.reference(u,T0)
    setup(max(T,T0),k,T,seed.eta,3,2,2,8,32.,.82,cache=panels,cache_only=True)
    model=Late(k,H,panels,True,compact=True,local_jac=True)
    A=moments(p,th)
    S0=model.moment(u)[1]
    data={}
    for order in ("be","mid","sdirk2"):
        y=[]
        for n in N:
            q=u.copy()
            e=0.;iterations=0;failed=None
            try:
                for _ in range(n):
                    h=step/n
                    if order=="be":
                        q,info=implicit(model,q,h,T)
                    elif order=="mid":
                        q,info=midpoint(model,q,h,T,p)
                    else:
                        q,info=sdirk2(model,q,h,T)
                    e=max(e,info["residual"]); iterations+=info["newton_iterations"]
            except (RuntimeError,FloatingPointError,RuntimeWarning,ValueError,AssertionError) as exc:
                failed=str(exc) or type(exc).__name__
            if failed:
                y.append(dict(n=n,rejected=failed)); continue
            B=moments(p,model.arrays(q))
            err=max(abs(B[:,0].sum(axis=1)/A[:,0].sum(axis=1)-1.))
            t=e_total(B)/e_total(A)-1.
            assert max(err,abs(t))<3e-9
            y.append(dict(n=n,end=B.tolist(),residual=e,newton_iterations=iterations,
                          entropy_change=model.moment(q)[1]/S0-1.,count_error=float(err),energy_error=t))
        data[order]=y
    metadata=dict()
    for key,v in data.items():
        assert len(v)==3 and all("end" in r for r in v),data
        f=np.array([r["end"] for r in v])
        defect=np.sum(abs(np.diff(f,axis=0)),axis=-1)/np.sum(abs(f[-1]),axis=-1)
        metadata[key]=dict(weighted_shell_L1_adjacent_defects=defect.tolist(),
                           observed_log2_ratio=np.log2(defect[0]/defect[1]).tolist())
    assert np.max(metadata["sdirk2"]["weighted_shell_L1_adjacent_defects"][-1])<6e-7
    assert np.min(metadata["mid"]["observed_log2_ratio"])<0.
    return dict(shape_sha256=sha(shape),eos_sha256=sha(eos),frozen_log_expansion_microstep_size=step,
        methods=data,metric=metadata,
        scope="Frozen 12 MeV collision-substep comparison; no FRW work, pair kicks, source, or full-trajectory convergence. Reflected midpoint is tested but not selected for this stiff problem. The two-stage solver uses the same moment fluxes; it is not unconditionally positivity-preserving.")

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos",type=Path)
    p.add_argument("shape",type=Path)
    p.add_argument("--step",type=float,default=.002)
    p.add_argument("--steps",type=int,nargs="+",default=[2,4,8])
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.shape,a.step,a.steps),indent=2,allow_nan=False))
