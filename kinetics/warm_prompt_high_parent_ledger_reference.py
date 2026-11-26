"""Measure high-parent on-shell prompt sources with explicit reacting support."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from joint_prompt_two_moment_evolution_reference import Collision
from joint_resonant_frw_split_reference import pack
from joint_visible_heavy_adiabatic_reference import angle_panels,snapshot_functions
from late_archived_shape_continuation_reference import Late,collide
from three_species_decay_onsager_reference import M,m,mp,width,P,Es
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom,pad


def high(s,cut,n):
    p=np.asarray(s["lab_momentum_edges_GeV"])
    phi,bp,chi,bc=snapshot_functions(s)
    threshold=np.hypot(M,cut)
    upper=bp[-1]+bc[-1]
    if upper<=threshold:
        return np.zeros(4)
    # Include the exact high-parent cut and a diagnostic on each
    # initial daughter at the half-cut. Both physical angular boundaries
    # and crossings of all daughter-cell intervals contribute panels.
    tag=cut/2
    ephi,echi=np.hypot(mp,tag),np.hypot(m,tag)
    sp=np.unique(np.r_[bp,ephi]);sc=np.unique(np.r_[bc,echi])
    def gate(E,d,E0):
        px=np.sqrt(np.maximum((E-d)*(E+d),0.))
        hi=E*E0+px*P
        return np.r_[M*(E0*E0+px*px)/hi,M*hi/(d*d)]
    edge=np.unique(np.clip(np.r_[threshold,upper,gate(sp,mp,Es),gate(sc,m,M-Es),
                                 (sp[:,None]+sc).ravel()],threshold,upper))
    z,w=roots_legendre(n)
    interval=np.diff(edge)/2
    E=(((edge[:-1]+edge[1:])/2)[:,None]+interval[:,None]*z).ravel()
    W=(interval[:,None]*w).ravel()
    pH=np.sqrt((E-M)*(E+M))
    a=np.zeros(4)
    for j in range(0,len(E),160):
        ii=slice(j,j+160)
        eps,w=angle_panels(E[ii],sp,sc,angle=n)
        el=E[ii,None,None]-eps
        f,g=phi(eps),chi(el)
        A=w*f*g
        B=w*(1+f)*(1-g)
        fH=A.sum(axis=(1,2))/(A+B).sum(axis=(1,2))
        geom=(2/np.pi**2)*pH[ii]*W[ii]
        live=geom*width*M*(1-fH)
        tagged_e=(A*(eps*(eps>=ephi)+el*(el>=echi))).sum(axis=(1,2))
        a+=[live@A.sum(axis=(1,2)),live@(A.sum(axis=(1,2))*E[ii]),
            live@tagged_e,geom@(E[ii]**2*fH)]
    assert a[2]>=a[1]/2*(1-4e-14)
    assert a[2]<=a[1]*(1+4e-14)
    return a


def run(D):
    a=read(D/"warm-light-occupation-next-q32-tape-2026-11-26.json.gz")[0]
    S=a["snapshot"]
    assert digest(S)==a["snapshot_sha256"]
    p=np.asarray(S["lab_momentum_edges_GeV"])
    q=np.r_[p,p[-1]*np.geomspace(1.,2.,6)[1:]]
    H=read(D/"warm-near-rest-light-inverse-decay-2026-11-26.json.gz")["cases"][1]["stages"][0]["fine"]["H_over_rest_decay_width"]*width
    th=pad(S,q); v=pack(th)
    rho=get_mom(q,th)[2].sum()
    model=Late(Collision(q),H,[],True,local_jac=True)
    rows=[]
    def add(t,A):
        ss=dict(lab_momentum_edges_GeV=q,entropy_variable_scalar=A[0],entropy_variable_fermion=A[1])
        vals=np.asarray([high(ss,1.,n) for n in (6,12,18)])
        y=vals.copy()
        y[:,:3]/=rho*H
        y[:,0]*=M  # use heavy rest energy to normalize count flux
        y[:,3]/=rho
        if t==0.:
            assert np.all(y==0.)
        else:
            assert np.all(y>0.) and np.max(abs(y[-1]/y[-2]-1))<.02
        rows.append(dict(elapsed_times_H=t,quadrature_6_12_18=y.tolist()))
    add(0.,th)
    times=[0.,1e-8,5e-8,2e-7,1e-6,5e-6]
    for a,b in zip(times[:-1],times[1:],strict=True):
        v,*_=collide(model,v,b-a,S["T_visible_GeV"],continuation=True,balanced=True)
        add(b,model.arrays(v))
    return dict(scope="Fixed-clock explicit q32-archive to vacuum-padded q64 prompt dynamics. "
                      "Only the bare H inverse-decay source, no other heavy production or exterior "
                      "heating. Snapshots not a cosmological upper or cumulative gross ledger.",
                columns=["high_inverse_number_times_M_over_rho_H","high_inverse_energy_over_rho_H",
                         "tagged_daughter_energy_in_high_inverse_over_rho_H","high_frozen_heavy_energy_over_retained_rho"],
                rows=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--tapes",type=Path,default=Path("results"))
    print(json.dumps(run(p.parse_args().tapes),indent=2,allow_nan=False))
