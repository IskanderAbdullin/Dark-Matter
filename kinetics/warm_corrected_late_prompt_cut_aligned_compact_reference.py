"""Close late prompt evolution on shells aligned to the excess-energy cut."""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from joint_prompt_two_moment_evolution_reference import Collision,kin
from joint_resonant_frw_split_reference import pack
from joint_resonant_centered_stiff_reference import advance
from late_archived_shape_continuation_reference import collide
from late_radiation_background_reference import MPL, TabulatedRadiation
from three_species_decay_onsager_reference import M,m,mp
from warm_corrected_exterior_reactivated_tail_stress_reference import reactivate
from warm_coupled_serial_restart_reference import atomic,sha
from warm_compact_absolute_prompt_reference import LocalAbsolutePrompt as DirectPrompt
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_excess_energy_tail_reference import integrate
from warm_prompt_reentry_shape_audit_reference import get_mom
from joint_visible_heavy_adiabatic_reference import snapshot_functions


def excess(state,cut):
    p=np.asarray(state["lab_momentum_edges_GeV"])
    B,_,F,_=snapshot_functions(state)
    z,w=roots_legendre(16)
    values=[]
    for d,g,f in ((mp,1,B),(m,4,F)):
        threshold=np.sqrt((cut-d)*(cut+d))
        e=np.unique(np.r_[max(threshold,p[0]),p[p>threshold]])
        q=e[:-1,None]+np.diff(e)[:,None]*(1+z)/2
        E=np.hypot(q,d)
        values.append(float(g/(4*np.pi**2)*np.sum(np.diff(e)[:,None]*w*q*q*f(E)*(E-cut))))
    return values


def split_at_cut(A,cut):
    p=np.asarray(A["lab_momentum_edges_GeV"])
    edges=[np.sqrt((cut-d)*(cut+d)) for d in (mp,m)]
    r=np.unique(np.r_[p,edges])
    k=np.searchsorted(p,(r[:-1]+r[1:])/2,side="right")-1
    assert np.all((k>=0)&(k<len(p)-1))
    assert np.all(r[:-1]>=p[k]) and np.all(r[1:]<=p[k+1])
    D=dict(A,lab_momentum_edges_GeV=r.tolist())
    for d,sp in ((mp,"scalar"),(m,"fermion")):
        x=kin(p,d);y=kin(r,d);w=np.diff(x)
        key="entropy_variable_"+sp
        a,b=np.asarray(A[key])
        D[key]=[ (a[k]+b[k]*(y[:-1]-x[k])/w[k]).tolist(), (b[k]*np.diff(y)/w[k]).tolist() ]
    return D


def observable_weights(model,C):
    pieces=[]
    for c in model.k.c:
        lo=c.low+c.mass
        hi=c.hi+c.mass
        assert np.all((lo>=C-2e-16)|(hi<=C+2e-16))
        take=hi>C+2e-16
        pieces.append(np.stack((np.where(take,np.maximum(lo-C,0.),0.),
                                np.where(take,c.dk,0.)),axis=1).ravel())
    return np.concatenate(pieces)


def predictor(model,x,h):
    target=model.moment(x)[0].reshape(2,model.n,2)
    y=x.copy()
    for _ in range(5):
        F=model.collision(y)[0].reshape(2,model.n,2)
        th=[a.copy() for a in model.arrays(y)]
        for k,c in enumerate(model.k.c):
            old=target[k]
            est=old+h*F[k]
            # First-order positive targets just for bins about to gain
            # significantly. Keep the physical implicit target intact.
            use=(est[:,0]>1.25*old[:,0])&(est[:,0]>1e-205)
            use &= (est[:,1]>0.)&(est[:,1]<est[:,0])
            if not np.any(use):
                continue
            N=np.where(use,est[:,0],0.)
            U=np.where(use,est[:,1],0.)
            try:
                a,_,_=c.reconstruct((N,c.low*N+c.dk*U))
            except (ValueError,AssertionError,RuntimeError,FloatingPointError):
                continue
            th[k][:,use]=np.asarray(a)[:,use]
        y=pack(th)
    return y


def run(eos,prior,archive,out,hard,n,seed_floor,guided):
    v=read(archive)
    S=v["snapshot"]; assert v["snapshot_sha256"]==digest(S)
    T=S["T_visible_GeV"]
    cut=np.hypot(M,1./read(prior)["branches"][-1]["entropy_expansion_from_original_anchor"])
    A,_,_,_,_=reactivate(S,hard,.005,cap=1.05)
    A=split_at_cut(A,cut)
    keys=("entropy_variable_scalar","entropy_variable_fermion")
    if seed_floor:
        assert 0.<seed_floor<1e-60
        A=dict(A)
        for k in keys:
            v=np.asarray(A[k])
            empty=(v[0]>=700.)&(v[1]==0.)
            v[0,empty]=-np.log(seed_floor)
            A[k]=v.tolist()
    p=np.asarray(A["lab_momentum_edges_GeV"])
    N,_,E=get_mom(p,[A[k] for k in keys])
    rho=E.sum()
    H=np.sqrt((TabulatedRadiation(eos,T).state(1.)[1]+rho)/3)/MPL
    model=DirectPrompt(Collision(p),H,[],True)
    weights=observable_weights(model,cut)
    plan=dict(eos_sha256=sha(eos),seed_sha256=digest(A),method_sha256=sha(Path(__file__)),compact_sha256=sha(Path(__file__).with_name("warm_compact_absolute_prompt_reference.py")),hard_GeV=hard,substeps=n,artificial_empty_shell_filling=seed_floor,positive_predictor=guided)
    out.mkdir(parents=True,exist_ok=True)
    rows=[]
    t=np.r_[0.,np.geomspace(1e-11,1e-6,6)]
    x=pack([A[k] for k in keys])
    for i,b in enumerate(t):
        path=out/f"step-{i}.json"
        if path.exists():
            R=read(path)
            assert R["plan"]==plan and R["time"]==b and R["hash"]==digest(R["snapshot"])
            x=pack([R["snapshot"][k] for k in keys])
        else:
            if i:
                for j in range(n):
                    h=(b-t[i-1])/n
                    if guided:
                        try:
                            guide=predictor(model,x,h)
                            x,*_=advance(model,x,h,T,initial_theta=guide,vacuum_scaled=True,balanced_merit=True)
                            continue
                        except (ValueError,AssertionError,RuntimeError,FloatingPointError):
                            pass
                    x,*_=collide(model,x,h,T,continuation=True,balanced=True)
            st=dict(A,**{k:a.tolist() for k,a in zip(keys,model.arrays(x),strict=True)})
            m,_,F=get_mom(p,[st[k] for k in keys])
            defect=[float(np.max(abs(m.sum(axis=1)/N.sum(axis=1)-1.))),float(abs(F.sum()/rho-1.))]
            assert max(defect)<6e-12
            native=weights@model.moment(x)[0]/rho
            measured=float(sum(excess(st,cut))/rho)
            assert abs(native-measured)<5e-11*max(native,measured,1e-60)
            R=dict(plan=plan,hash=digest(st),snapshot=st,time=b,moment_defects=defect,
                   native_signed_cut_flux=float(weights@model.collision(x)[0]/rho),
                   cut_moment_error_over_rho=float(abs(native-measured)),
                   excess_over_initial_energy=(np.asarray(excess(st,cut))/rho).tolist(),
                   independent_closed_prompt_gross_and_signed=(integrate(st,cut,8,True)/(rho*H)).tolist())
            atomic(path,R)
        rows.append({k:a for k,a in R.items() if k!="snapshot"})
        print(json.dumps(dict(done=i,hard=hard,substeps=n)),file=sys.stderr,flush=True)
    return dict(scope="Fixed-clock prompt-only stress on finite closed padded support, not late FRW or untracked shape.",
                rows=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","prior","archive","out"):
        p.add_argument(k,type=Path)
    p.add_argument("--hard",type=float,required=True)
    p.add_argument("--substeps",type=int,required=True)
    p.add_argument("--seed-floor",type=float,default=0.)
    p.add_argument("--predictor",action="store_true")
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.prior,a.archive,a.out,a.hard,a.substeps,a.seed_floor,a.predictor),indent=2,allow_nan=False))
