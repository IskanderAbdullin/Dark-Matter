"""Scalar pair gain and loss at exactly zero tagged light momentum."""
import argparse
import json
import numpy as np
from pathlib import Path
from scipy.special import expit,roots_legendre
from archived_scalar_pair_reaction_reference import amplitude
from joint_visible_heavy_adiabatic_reference import snapshot_functions
from warm_light_occupation_stage_audit_reference import digest,read
from three_species_decay_onsager_reference import m,mp,width

def calc(B,be,F,fe,r,a,thermal=False):
    Emax=max(fe[-1],2*be[-1]-m)
    edge=np.unique(np.clip(np.r_[m,fe,2*be-m,Emax],m,Emax))
    x,w=roots_legendre(r); z,zw=roots_legendre(a)
    E=(((edge[:-1]+edge[1:])/2)[:,None]+np.diff(edge)[:,None]*x/2).ravel()
    dE=(np.diff(edge)[:,None]*w/2).ravel()
    p=np.sqrt((E-m)*(E+m))
    N=p*E*dE/np.pi**2
    s=2*m*(m+E)
    b=np.sqrt(s/4-mp*mp)
    center=(m+E)/2
    D=p*b/np.sqrt(s)
    # Polar cosine wrt boost, initial tagged L is antiparallel.
    gate=(be-center[:,None])/D[:,None]
    cut=np.sort(np.clip(np.concatenate([-np.ones((len(E),1)),gate,-gate,np.ones((len(E),1))],axis=1),-1,1),axis=1)
    d=np.diff(cut,axis=1)/2
    pol=(((cut[:,1:]+cut[:,:-1])/2)[:,:,None]+d[:,:,None]*z).reshape(len(E),-1)
    wp=(d[:,:,None]*zw).reshape(len(E),-1)
    f0=F(np.array([m]))[0]; f1=F(E)[:,None]
    b0=B(center[:,None]+D[:,None]*pol)
    b1=B(center[:,None]-D[:,None]*pol)
    G=(1-f0)*(1-f1)*b0*b1
    L=f0*f1*(1+b0)*(1+b1)
    Q=(N*b/(64*np.pi*m*E*np.sqrt(s)))[:,None]*wp*amplitude(s[:,None],-pol)
    if thermal:
        defect=np.max(abs(G-L)/np.maximum(np.maximum(G,L),1e-280))
        assert defect<5e-14
    else:
        defect=None
    return np.array([np.sum(Q*G),np.sum(Q*L)]),defect


def verify():
    T,mu=.022,-.02
    B=lambda E:1/np.expm1((E-mu)/T)
    F=lambda E:expit(-(E-mu)/T)
    be=np.hypot(np.linspace(0.,1.5,21),mp)
    fe=np.hypot(np.linspace(0.,1.5,21),m)
    x,d=calc(B,be,F,fe,10,9,thermal=True)
    defect=abs(x[0]/x[1]-1)
    assert defect<1e-14
    return dict(thermal_relative_integrand_defect=d,thermal_integrated_defect=defect)


def run(D):
    ref=read(D/"warm-near-rest-light-inverse-decay-2026-11-26.json.gz")
    out=[]
    q16=[]
    for q in (16,32):
        A=read(D/f"warm-light-occupation-next-q{q}-tape-2026-11-26.json.gz")
        A+=[z["checkpoint"] for z in read(D/f"warm-observed-suboperator-q{q}-tape-2026-11-26.json.gz")["blocks"]]
        matched=next(z for z in ref["cases"] if z["q"]==q)["stages"]
        rows=[]
        for j,(r,S0) in enumerate(zip(A,matched,strict=True)):
            S=r["snapshot"]
            assert r["snapshot_sha256"]==S0["snapshot_sha256"]==digest(S)
            B,be,F,fe=snapshot_functions(S)
            H=S0["fine"]["H_over_rest_decay_width"]*width
            stages=[]
            for order in ((8,8),(14,12),(24,16)):
                x,_=calc(B,be,F,fe,*order)
                x=x/H
                stages.append(dict(order=order,gain=x[0],loss=x[1],net=x[0]-x[1]))
            assert stages[-1]["net"]<0
            assert abs(stages[-1]["net"]/stages[-2]["net"]-1)<3e-5
            row=dict(T_GeV=S0["T_GeV"],snapshot_sha256=r["snapshot_sha256"],controls=stages)
            if q==16:
                q16.append(be[-1])
            else:
                Emax=q16[j]
                cropped=lambda E:np.where(E<Emax,B(E),0.)
                v,_=calc(cropped,np.r_[be[be<Emax],Emax],F,fe,24,16)
                row["same_q32_only_scalar_crop"]=[float(v[0]/H),float(v[1]/H),float((v[0]-v[1])/H)]
            rows.append(row)
        out.append(dict(q=q,stages=rows))
    return dict(scope="Frozen retained two-body scalar pair channel at exactly rest for the tagged light "
                      "spin and charge. No integrated relic or evolving number inference.",
                controls=verify(),cases=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tapes",type=Path,default=Path("results"))
    print(json.dumps(run(p.parse_args().tapes),indent=2,allow_nan=False))
