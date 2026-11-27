"""Matched stationary and finite-parent prompt evolution on a comoving grid."""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

from frw_exterior_cohort_reference import Exterior
from joint_resonant_centered_stiff_reference import advance
from joint_resonant_frw_split_reference import pack,redshift
from late_radiation_background_reference import MPL,TabulatedRadiation
from three_species_decay_onsager_reference import M
from warm_compact_absolute_prompt_reference import LocalAbsolutePrompt
from warm_comoving_parent_chart_reference import ComovingCollision,common_rule,remapless_parent_work
from warm_comoving_parent_competing_reference import FiniteParentElastic,StationaryParentElastic,midpoint_pair_kick
from warm_corrected_late_prompt_cut_aligned_reference import excess,predictor
from warm_coupled_serial_restart_reference import atomic,sha
from warm_late_finite_parent_prompt_reference import FiniteParent
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom,pad

KEYS=("entropy_variable_scalar","entropy_variable_fermion")


def parent_equilibrium(k,th):
    f=k.c[0].f(th[0],k.ph,k.I)
    g=k.c[1].f(th[1],k.ch,k.J)
    a=np.bincount(k.y,weights=k.pre*f*g,minlength=k.n)
    b=np.bincount(k.y,weights=k.pre*(1.+f)*(1.-g),minlength=k.n)
    return np.divide(a,a+b,out=np.zeros_like(a),where=a+b>0.)


def run(eos,paths,out,mode,order,common_subdivisions,outer=0.,restrict=False,root_tol=8e-12):
    R=[read(r) for r in paths]
    assert all(x["hash"]==digest(x["snapshot"]) and x["plan"]==R[0]["plan"] for x in R)
    assert R[0]["plan"]["eos_sha256"]==sha(eos)
    S=R[0]["snapshot"]
    T=S["T_visible_GeV"]
    q=np.asarray(S["lab_momentum_edges_GeV"])
    assert all(np.array_equal(r["snapshot"]["lab_momentum_edges_GeV"],q) for r in R)
    assert all(r["snapshot"]["T_visible_GeV"]==T for r in R)
    old_shells=len(q)-1
    if outer:
        assert outer>q[-1] and mode in ("expand-finite","expand-stationary")
        q=np.r_[q,np.geomspace(q[-1],outer,4)[1:]]
    else:
        assert not restrict
    C=R[0]["plan"]["normalized_cut_GeV"]
    times=np.array([r["time"]-R[0]["time"] for r in R])
    assert np.all(np.diff(times)>0.)
    zero=(mode.startswith("frozen-"))
    actual=np.exp(times) if not zero else np.ones(len(times))
    quad_scales=np.exp(np.linspace(0.,times[-1],2*common_subdivisions+1))
    Q,dN=common_rule(q,quad_scales,order)
    th=pad(S,q) if outer else [np.asarray(S[k]) for k in KEYS]
    op0=ComovingCollision(q,1.,Q,dN,retained_shells=old_shells if outer else None)
    h0=parent_equilibrium(op0,th)
    hv=h0.copy()
    del op0
    fin="finite" in mode
    extra=mode.endswith("-all")
    ex=Exterior()
    outside,exchange=np.zeros(2),np.zeros(2)
    max_pair_defect=0.
    assert not S["positive_exterior_restart"]["cohorts"]
    ta=TabulatedRadiation(eos,T)
    work=0.
    n,_,e=get_mom(q,th)
    rho=e.sum()
    initial_N=n.sum(axis=1)+(dN@hv if fin else 0.)
    initial_E=rho+((dN*np.hypot(M,Q))@hv if fin else 0.)
    out.mkdir(exist_ok=True,parents=True)
    rows=[]
    for j in range(1,len(R)):
        a,c=actual[j-1:j+1]
        b=np.sqrt(a*c)
        h=times[j]-times[j-1]
        for_stage=ta.state(b)
        def drift(a,b,theta):
            if a==b:
                return theta,0.
            u,w,*_=redshift(q,a,b,theta)
            return u,float(w.sum())
        th,dw=drift(a,b,th)
        work+=dw
        if fin:
            work+=remapless_parent_work(Q,dN,hv,a,b)["pressure_trace_work"]
        p=q/b
        # Hold the gravitating parent equal to the same collisionless
        # initial phase in both branches. The localized reaction changes
        # its energy ledger but not this clock prescription.
        grav=(dN*np.hypot(M,Q/b))@h0/b**3
        rho_mid=get_mom(p,th)[2].sum()
        ex.advance(b)
        H=np.sqrt((for_stage[1]+rho_mid+grav+ex.energy(b)/b**4)/3.)/MPL
        if extra:
            th,change,new,err=midpoint_pair_kick(p,th,b,h/(2*H),ex)
            outside+=new
            exchange+=change
            max_pair_defect=max(max_pair_defect,err)
        k=ComovingCollision(q,b,Q,dN,retained_shells=old_shells if restrict else None)
        if fin:
            op=FiniteParentElastic(k,H,for_stage[0]) if extra else FiniteParent(k,H)
            np.testing.assert_allclose(op.measure,dN/b**3,rtol=2e-14,atol=1e-300)
            op.old=hv
            op.dt=h
            ent=lambda u,after,before: after+op.heavy_entropy(op.occupation(u))>=before+op.heavy_entropy(hv)
        else:
            op=StationaryParentElastic(k,H,for_stage[0]) if extra else LocalAbsolutePrompt(k,H,[],True,stream_local=True)
            ent=None
        target=pack(th)
        gs=[pack(drift(1.,b,pad(r["snapshot"],q) if outer else
                       [np.asarray(r["snapshot"][key]) for key in KEYS])[0]) for r in R[j-1:j+1]]
        guess=target+(gs[1]-gs[0])
        if outer and not restrict:
            seeded=predictor(op,target,h)
            ph=op.arrays(guess)
            # Only guesses, and only in the additional cells. The
            # predecessor itself retains exact zero occupations there.
            for g,s in zip(ph,op.arrays(seeded),strict=True):
                g[:,old_shells:]=s[:,old_shells:]
            guess=pack(ph)
        def log(i,r,index):
            if i in (0,1,2,5,10,20,40,80):
                print(mode, j, i, f"{r:.5g}", index, file=sys.stderr,flush=True)
        v,stat=advance(op,target,h,for_stage[0],initial_theta=guess,vacuum_scaled=True,
                      balanced_merit=True,entropy_guard=ent,observer=log,residual_tolerance=root_tol)
        if fin:
            newparent=op.occupation(v)
            assert np.all((0.<=newparent)&(newparent<=1.))
            hv=newparent
        th=op.arrays(v)
        if extra:
            th,change,new,err=midpoint_pair_kick(p,th,b,h/(2*H),ex)
            outside+=new
            exchange+=change
            max_pair_defect=max(max_pair_defect,err)
        ex.advance(c)
        th,dw=drift(b,c,th)
        work+=dw
        if fin:
            work+=remapless_parent_work(Q,dN,hv,b,c)["pressure_trace_work"]
        n,_,e=get_mom(q/c,th)
        N=n.sum(axis=1)*c**3+(dN@hv if fin else 0.)
        E=c**4*e.sum()+(c*((dN*np.hypot(M,Q/c))@hv) if fin else 0.)
        R_N=(N+outside-exchange)/initial_N-1.
        R_E=(E+ex.energy(c)-work-ex.trace)/initial_E-1.
        assert abs(sum(exchange))/sum(initial_N)<3e-12
        assert max(max(abs(R_N)),abs(R_E))<3e-10
        P=dict(S,lab_momentum_edges_GeV=(q/c).tolist(),scale_factor=S["scale_factor"]*c,
                T_visible_GeV=ta.state(c)[0],
                represented_species_count_a3=(n.sum(axis=1)*(S["scale_factor"]*c)**3).tolist(),
                **{key:np.asarray(z).tolist() for key,z in zip(KEYS,th,strict=True)})
        def d(y):
            return hashlib.sha256(np.asarray(y,dtype="<f8").tobytes()).hexdigest()
        r=dict(a=float(a),b=float(b),c=float(c),time_J=R[j]["time"],step_H=float(h),H_GeV=float(H),
               initial_daughter_hash=R[0]["hash"],daughter_hash=digest(P),parent_hash=d(hv) if fin else None,
               W_over_initial_daughter_rho=float(sum(excess(P,C))/rho),solver=stat,
               pressure_trace_work_over_energy=float(work/initial_E),relative_N_defects=R_N.tolist(),
               relative_energy_defect=float(R_E),qsha256=d(Q),wsha256=d(dN),parent_nodes=len(Q),
               active_angular_nodes=len(k.y))
        if outer:
            from three_species_decay_onsager_reference import m,mp
            # New cells are wholly above the physical cut in these tests.
            assert min(np.hypot(q[old_shells]/c,x) for x in (mp,m))>=C
            newN=n[:,old_shells:].sum(axis=1)
            newE=e[:,old_shells:].sum(axis=1)
            if restrict:
                assert np.all(newN==0.) and np.all(newE==0.)
            r.update(daughter_original_outer_edge=p[old_shells],
                outer_prompt_is_open=not restrict,
                newly_resolved_excess_by_species_over_rho=((newE-C*newN)/rho).tolist(),
                newly_resolved_energy_by_species_over_rho=(newE/rho).tolist(),
                newly_resolved_number_by_species_a3=(newN*c**3).tolist())
        if extra:
            r.update(pair_comoving_count_exchange_over_initial=(exchange/initial_N).tolist(),
                     pair_outside_energy_over_initial=float(ex.energy(c)/initial_E),
                     pair_conservation_defect=max_pair_defect,elastic_panels=len(op.panels))
        rows.append(r)
        atomic(out/f"step-{j}.json",dict(row=r,daughter_snapshot=P))
        del k,op
    start=float(sum(excess(S,C))/rho)
    scope="Source-free closed prompt channel on prescribed FRW clock; finite parent gravity uses matched passive prescription."
    if extra:
        scope+=" Includes selected closed subpole elastic and open-outgoing pair channels. Exterior does not reenter."
    atomic(out/"summary.json",dict(scope=scope,
                                  input_sha256=R[0]["hash"],eos_sha256=sha(eos),mode=mode,radial_order=order,
                                  cut_GeV=C,W_initial_over_rho=start,rows=rows))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path);p.add_argument("out",type=Path)
    p.add_argument("mode",choices=("expand-finite","expand-stationary","frozen-finite","frozen-stationary","expand-finite-all","expand-stationary-all"))
    p.add_argument("paths",type=Path,nargs="+")
    p.add_argument("--order",type=int,default=3)
    p.add_argument("--common-subdivisions",type=int,default=1)
    p.add_argument("--prompt-pad",type=float,default=0.)
    p.add_argument("--restrict-prompt-pad",action="store_true")
    p.add_argument("--root-tol",type=float,default=8e-12)
    a=p.parse_args()
    run(a.eos,a.paths,a.out,a.mode,a.order,a.common_subdivisions,a.prompt_pad,a.restrict_prompt_pad,a.root_tol)
