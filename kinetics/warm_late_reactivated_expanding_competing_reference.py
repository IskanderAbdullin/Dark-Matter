"""One matched late expansion block with optional elastic and signed pair reactions."""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from compact_elastic_panel_reference import compile_panels, compact_elastic
from frw_exterior_restart_reference import export_restart, import_restart
from joint_below_pole_elastic_linear_modes_reference import setup
from joint_prompt_two_moment_evolution_reference import Collision
from joint_resonant_centered_stiff_reference import advance
from joint_resonant_frw_split_reference import pack, redshift
from joint_visible_heavy_adiabatic_reference import snapshot_functions
from late_radiation_background_reference import MPL, TabulatedRadiation
from scalar_pair_expanding_continuation_reference import MASS, e_total, moments, reconstruct
from warm_compact_absolute_prompt_reference import LocalAbsolutePrompt
from warm_corrected_late_prompt_cut_aligned_reference import excess
from warm_coupled_serial_restart_reference import atomic, sha
from warm_light_occupation_stage_audit_reference import digest, read
from warm_pair_common_fermion_chart_reference import rate


KEYS=("entropy_variable_scalar","entropy_variable_fermion")


class PromptElastic(LocalAbsolutePrompt):
    """Use the absolute prompt flux and the closed compact elastic derivative."""
    def __init__(self, kernel, H, panels):
        super().__init__(kernel,H,[],True)
        self.panels=panels
        self.compiled=compile_panels(panels,self.n)

    def centered(self,z,ref,R0):
        a,da=self.collision(ref+z)
        if not self.panels:
            return a,da
        b,db=compact_elastic(z,self.t,self.eta,self.H,self.compiled,self.d,True)
        return a+b,lambda v: da(v)+db(v)


def phase(s, p, theta, b, T):
    n=moments(p,theta)[:,0].sum(axis=1)
    return dict(s,lab_momentum_edges_GeV=p.tolist(),**{
        k:np.asarray(v).tolist() for k,v in zip(KEYS,theta,strict=True)},
                scale_factor=s["scale_factor"]*b,T_visible_GeV=T,
                represented_species_count_a3=(n*(s["scale_factor"]*b)**3).tolist())


def run(eos,start,guide,mode,extent,predecessor=None):
    assert mode in ("prompt","prompt-pair","all")
    a,b=read(start),read(guide)
    A,G=a["snapshot"],b["snapshot"]
    assert all(z["hash"]==digest(z["snapshot"]) for z in (a,b))
    assert b["plan"]==a["plan"]
    assert a["plan"]["eos_sha256"]==sha(eos)
    raw_p=np.asarray(A["lab_momentum_edges_GeV"])
    assert np.array_equal(raw_p,G["lab_momentum_edges_GeV"])
    S=A if predecessor is None else predecessor["snapshot"]
    if predecessor is not None:
        assert predecessor["hash"]==digest(S) and predecessor["raw_target_hash"]==a["hash"]
        assert predecessor["mode"]==mode and predecessor["eos_sha256"]==sha(eos)
        assert np.allclose(raw_p*(A["scale_factor"]/S["scale_factor"]),
                           S["lab_momentum_edges_GeV"],rtol=4e-15,atol=0.)
    t=S["T_visible_GeV"]
    C=a["plan"]["normalized_cut_GeV"]
    q=np.array(S["lab_momentum_edges_GeV"])
    h=b["time"]-a["time"]
    beta,gamma=np.exp([h/2,h])
    p=q/beta
    clock=TabulatedRadiation(eos,t)
    T=clock.state(beta)[0]
    th=[np.asarray(S[k]) for k in KEYS]
    initial=moments(q,th)
    E0=e_total(initial)
    N0=initial[:,0].sum(axis=1)
    th,w1,*_=redshift(q,1.,beta,th)
    rho=e_total(moments(p,th))
    x0=import_restart(S["positive_exterior_restart"])
    e_ex0,trace0=x0.energy(1.),x0.trace
    x0.advance(beta)
    H=np.sqrt((clock.state(beta)[1]+rho+x0.energy(beta)/beta**4)/3.)/MPL
    max_pair_defect=0.
    exn=np.zeros(2)
    dn=np.zeros(2)
    pair_out_above=0.
    stages=[]
    def record(label,p0,th0,step):
        st=phase(S,p0,th0,step,clock.state(step)[0])
        W=sum(excess(st,C))/E0
        stages.append(dict(stage=label,relative_expansion=float(step),W_over_initial_rho=float(W),
                           sha=digest(st)))
    def kick(theta):
        nonlocal max_pair_defect,pair_out_above
        elapsed=h/(2*H)
        current=moments(p,theta)
        def rhs(theta,log_out=False):
            B,_,F,_=snapshot_functions(phase(S,p,theta,beta,T))
            result={}
            net,_,_,err,_=rate(p,F,B,(3,7,6,2),**(
                dict(excess_cut=(C,),excess_result=result) if log_out else {}))
            return net,err,result
        first,_,_=rhs(theta)
        mid=reconstruct(p,current+elapsed*.5*first[:,:,1:-1])
        second,err,W=rhs(mid,True)
        max_pair_defect=max(max_pair_defect,err)
        pair_out_above+=elapsed*W["rate_GeV5"][0][3]/E0
        increment=elapsed*second
        dn[:]+=beta**3*increment[:,0,1:-1].sum(axis=1)
        for j,mi in enumerate(MASS):
            for side,ix in enumerate((0,-1)):
                amount=beta**3*increment[j,0,ix]
                exn[j]+=amount
                x0.add(beta,mi,j,side,p[ix],amount,beta**4*increment[j,1,ix])
        v=reconstruct(p,current+increment[:,:,1:-1])
        return v
    record("first_redshift",p,th,beta)
    if mode != "prompt":
        th=kick(th)
        record("first_pair_half",p,th,beta)
    kernel=Collision(p)
    panels=[]
    if mode=="all":
        setup(t,kernel,T,(0.,0.),3,2,2,8,extent,.82,cache=panels,cache_only=True)
    m=PromptElastic(kernel,H,panels)
    target=pack(th)
    guess,_w,*_=redshift(q,1.,beta,[np.array(G[k]) for k in KEYS])
    # Transfer the exact prompt-only update in entropy coordinates as a
    # guess on the actual predecessor. The implicit target is unchanged.
    start_half,_w,*_=redshift(q,1.,beta,[np.array(A[k]) for k in KEYS])
    moved_guess=target+pack(guess)-pack(start_half)
    def log(i,r,j):
        if i in (0,1,2,5,10,20,40,80):
            print(f"{mode} {i=} {r=:.4g} {j=}",file=sys.stderr,flush=True)
    v,solver=advance(m,target,h,T,initial_theta=moved_guess,
                     vacuum_scaled=True,balanced_merit=True,observer=log)
    th=m.arrays(v)
    record("accepted_local_collision",p,th,beta)
    if mode != "prompt":
        th=kick(th)
        record("second_pair_half",p,th,beta)
    x0.advance(gamma)
    th,w2,*_=redshift(q,beta,gamma,th)
    st=phase(S,q/gamma,th,gamma,clock.state(gamma)[0])
    st["positive_exterior_restart"]=export_restart(x0)
    record("endpoint",q/gamma,th,gamma)
    final=moments(q/gamma,th)
    n=gamma**3*final[:,0].sum(axis=1)
    rel_count=(n-N0-dn)/N0
    overall_count=(n+exn).sum()/N0.sum()-1.
    energy=(gamma**4*e_total(final)+x0.energy(gamma)-e_ex0-w1.sum()-w2.sum()-(x0.trace-trace0))/E0-1.
    assert max(abs(overall_count),abs(energy),max(abs(rel_count)))<3e-10
    stages[-1]["sha"]=digest(st)
    return dict(input_hash=digest(S),raw_target_hash=b["hash"],guide_hash=b["hash"],snapshot=st,hash=digest(st),
                mode=mode,eos_sha256=sha(eos),relative_log_expansion=float(h),
                fixed_energy_cut_GeV=C,stages=stages,relative_energy_residual=float(energy),
                relative_total_count_residual=float(overall_count),
                per_species_relative_count_residual=rel_count.tolist(),
                pair_count_increment_over_start=(dn/N0).tolist(),pair_outgoing_W_at_collision_per_initial_rho=float(pair_out_above),
                max_pair_defect=max_pair_defect,number_elastic_panels=len(panels),
                elastic_extent_in_T=extent if panels else None,solver=solver)


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","start","guide","out"):
        p.add_argument(k,type=Path)
    p.add_argument("mode",choices=("prompt","prompt-pair","all"))
    p.add_argument("--extent",type=float,default=32.)
    p.add_argument("--predecessor",type=Path)
    a=p.parse_args()
    pred=read(a.predecessor) if a.predecessor else None
    atomic(a.out,run(a.eos,a.start,a.guide,a.mode,a.extent,pred))
