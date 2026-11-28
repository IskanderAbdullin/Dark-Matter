"""Audit the fixed-comoving parent geometry against the native stationary chart."""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from joint_prompt_two_moment_evolution_reference import Collision
from joint_resonant_frw_split_reference import redshift
from late_radiation_background_reference import MPL,TabulatedRadiation
from three_species_decay_onsager_reference import M
from warm_comoving_parent_chart_reference import ComovingCollision,common_rule,remapless_parent_work
from warm_coupled_serial_restart_reference import atomic,sha
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom


def phase(k,th,C):
    f=k.c[0].f(np.asarray(th[0]),k.ph,k.I)
    g=k.c[1].f(np.asarray(th[1]),k.ch,k.J)
    a=k.pre*f*g
    b=k.pre*(1.+f)*(1.-g)
    A=np.bincount(k.y,weights=a,minlength=k.n)
    B=np.bincount(k.y,weights=b,minlength=k.n)
    nt=(A+B)[k.y]
    fparent=np.divide(A,A+B,out=np.zeros_like(A),where=A+B>0.)
    numerator=b*A[k.y]-a*B[k.y]
    W=np.maximum(k.ph-C,0.)+np.maximum(k.ch-C,0.)
    flux=k.fac*np.divide(numerator,nt,out=np.zeros_like(nt),where=nt>0.)
    derived=np.array([flux.sum(), flux@k.ph, flux@k.ch, flux@W])
    gross=k.fac@(a*(1.-fparent[k.y])+b*fparent[k.y])
    assert abs(derived[0])<5e-13*gross
    assert abs(derived[1]+derived[2])<5e-13*gross
    return fparent,derived,gross


def run(eos,path,order,span):
    a=read(path)
    S=a["snapshot"]
    assert a["hash"]==digest(S) and a["plan"]["eos_sha256"]==sha(eos)
    q=np.array(S["lab_momentum_edges_GeV"])
    C=a["plan"]["normalized_cut_GeV"]
    T=S["T_visible_GeV"]
    theta=[np.asarray(S["entropy_variable_"+k]) for k in ("scalar","fermion")]
    rho=get_mom(q,theta)[2].sum()
    clock=TabulatedRadiation(eos,T)
    scale=np.exp(np.array([0.,.5,1.])*span)
    Q,dN=common_rule(q,scale,order)
    rows=[]
    f0=None
    for b in scale:
        th=theta if b==1. else redshift(q,1.,b,theta)[0]
        lab=q/b
        den=get_mom(lab,th)[2].sum()
        H=np.sqrt((clock.state(b)[1]+den)/3.)/MPL
        ka=Collision(lab)
        kb=ComovingCollision(q,b,Q,dN)
        xf,fx,mg=phase(ka,th,C)
        y,fy,mh=phase(kb,th,C)
        # Bounded positive occupations and background conservation are
        # separate from convergence of the detailed collision law.
        assert np.all((0.<=y)&(y<=1.))
        if f0 is None:
            f0=y
        Mtot=dN@(np.hypot(M,Q/b)*y)/b**3
        rows.append(dict(scale=float(b),H_GeV=float(H),parent_nodes=[ka.n,kb.n],
                         angular_samples=[len(ka.y),len(kb.y)],
                         normalized_gross_number=[mg/(rho*H),mh/(rho*H)],
                         stationary_parent_dark_energy_fraction=float(Mtot/rho),
                         native_and_common_cut_flux=[fx[3]/(rho*H),fy[3]/(rho*H)],
                         relative_common_vs_native_W_flux=float(fy[3]/fx[3]-1.),
                         common_number_and_energy_relative_residual=max(abs(fy[0]),abs(fy[1]+fy[2]))/mh))
        del ka,kb
        print(order,float(b),rows[-1]["relative_common_vs_native_W_flux"],file=sys.stderr,flush=True)
    work=remapless_parent_work(Q,dN,f0,scale[0],scale[-1])
    return dict(hard_GeV=a["plan"]["hard_GeV"],scope="Stationary geometry parity and pure parent Liouville transport; no finite reaction feedback.",
                phase_hash=a["hash"],reference_T_GeV=T,cut_GeV=C,order=order,span_log_a=span,rows=rows,
                collisionless_parent=work)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","input","output"):
        p.add_argument(k,type=Path)
    p.add_argument("--order",type=int,default=4)
    p.add_argument("--span",type=float,default=9.e-7)
    a=p.parse_args()
    atomic(a.output,run(a.eos,a.input,a.order,a.span))
