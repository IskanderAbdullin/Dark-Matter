"""Frozen tagged light-scalar elastic rates below a fixed invariant cut."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import expit, roots_legendre

from joint_visible_heavy_adiabatic_reference import snapshot_functions
from scalar_rest_yukawa_elastic_ir_reference import cross,trace_control,dot
from three_species_decay_onsager_reference import M,m,mp,y,width
from warm_light_occupation_stage_audit_reference import digest,read


def trace(s,xc):
    """Elastic differential cross section in d cosine / 2, vectorized."""
    s=np.asarray(s)
    V=np.sqrt(s)
    q=np.sqrt(np.maximum((s-(m+mp)**2)*(s-(m-mp)**2),0.))/(2*V)
    ef=(s+m*m-mp*mp)/(2*V)
    ep=V-ef
    z=np.zeros_like(xc)
    qx=q*np.sqrt(np.maximum(0.,1-xc*xc))
    p1=np.stack([ef+z,z,z,q+z],axis=-1)
    k2=np.stack([ep+z,qx,z,q*xc],axis=-1)
    P=np.stack([V+z,z,z,z],axis=-1)
    p2=P-k2
    U=p1-k2
    a=s-M*M
    b=dot(U,U)-M*M
    f=M*(1/a+1/b)
    Q=P/a[...,None]+U/b[...,None]
    q0=dot(p1,p2)
    A=2*f*f*(q0+m*m)+4*f*m*dot(p1+p2,Q)+2*(2*dot(p1,Q)*dot(p2,Q)+(m*m-q0)*dot(Q,Q))
    assert np.min(A)>=-1e-12
    return y**4*np.maximum(A,0)/(16*np.pi*s)


def collision(p,B,be,F,fe,eta,orders,details=False):
    e0=np.hypot(m,p)
    # For fixed e the smallest possible s obeys
    # s >= m^2 + mp^2 + 2*(e0-p)*e. This covers even a
    # partner outside the retained grid in the *actual final* state.
    Smax=(eta*M)**2
    L=(Smax-m*m-mp*mp)/2
    ecut=L*(e0+p)/(m*m)
    edges=np.unique(np.clip(np.r_[mp,be,ecut],mp,ecut))
    nodes,ww=roots_legendre(orders[0])
    mid=(edges[:-1]+edges[1:])/2
    de=np.diff(edges)/2
    incoming=(mid[:,None]+de[:,None]*nodes).ravel()
    incoming_w=(de[:,None]*ww).ravel()
    a,av=roots_legendre(orders[1])
    x,xw=roots_legendre(orders[2])
    cosaz=np.cos(np.pi*(np.arange(orders[3])+.5)/orders[3])
    q=np.sqrt((incoming-mp)*(incoming+mp))
    umin=np.clip(np.divide(e0*incoming-L,p*q,out=np.where(e0*incoming<=L,-1.,1.),
                           where=p*q>0.),-1.,1.)
    f0=F(np.array([e0]))[0]
    out=np.zeros(2)
    error=0.
    for j in np.flatnonzero(umin<1.):
        E=incoming[j]
        qq=q[j]
        wi=(1-umin[j])*av/4
        u=(1+umin[j])/2+(1-umin[j])*a/2
        S=m*m+mp*mp+2*(e0*E-p*qq*u)
        total=e0+E
        assert np.max(S)<=Smax*(1+3e-14)
        v=np.sqrt(S)
        kh=np.sqrt((S-(m+mp)**2)*(S-(m-mp)**2))/(2*v)
        Bst=np.sqrt(np.maximum(total*total-S,0.))
        B0=total*(S+mp*mp-m*m)/(2*S)
        W=Bst*kh/v
        midL=total-B0
        # Direction of the initial light leg relative to the total lab boost.
        if p==0.:
            polar=np.full_like(S,-1.)
        else:
            polar=np.clip((e0-midL)/np.maximum(W,1e-280),-1.,1.)
        # Outgoing scalar polar cosine measured relative to the lab boost.
        gate=np.concatenate([-np.ones((len(S),1)),(be-B0[:,None])/W[:,None],
                             (midL[:,None]-fe)/W[:,None],np.ones((len(S),1))],axis=1)
        gate=np.sort(np.clip(gate,-1.,1.),axis=1)
        d=np.diff(gate,axis=1)/2
        active=np.any(d>0,axis=0)
        center=((gate[:,1:]+gate[:,:-1])/2)[:,active]
        d=d[:,active]
        z=(center[:,:,None]+d[:,:,None]*x).reshape(len(S),-1)
        wa=(d[:,:,None]*xw/2).reshape(len(S),-1)
        mu=np.clip(polar[:,None,None]*z[:,:,None]+np.sqrt(np.maximum(0.,1-polar**2))[:,None,None]*
                    np.sqrt(np.maximum(0.,1-z*z))[:,:,None]*cosaz,-1.,1.)
        cross_avg=trace(S[:,None,None],mu).mean(axis=2)
        ep=B0[:,None]+W[:,None]*z
        ef=total-ep
        f1=F(ef);b1=B(ep);b0=B(np.array([E]))[0]
        G=(1-f0)*(1+b0)*f1*b1
        Ls=f0*b0*(1-f1)*(1+b1)
        if details:
            error=max(error,float(np.max(abs(G-Ls)/np.maximum(np.maximum(G,Ls),1e-280))))
        # p_phi^2 dp_phi = p_phi E_phi dE_phi, exact invariant
        # relative velocity sqrt((p0.p1)^2-m^2 mp^2)/(e0 E_phi).
        relative=np.sqrt(np.maximum(0.,(e0*E-p*qq*u)**2-(m*mp)**2))/(e0*E)
        measure=incoming_w[j]*qq*E/(2*np.pi**2)*wi*relative
        weight=wa*cross_avg*measure[:,None]
        out+=np.array([np.sum(weight*G),np.sum(weight*Ls)])
    return out,error


def controls():
    theta=np.linspace(-1.,1.,31)
    tests=[]
    for E in (.005,.01,.03,.07,.14,.19):
        s=m*m+mp*mp+2*m*E
        x=trace(s,theta)
        z=cross(0,E*m/mp,theta)
        tests.append(float(np.max(abs(x/z-1))))
    assert max(tests)<3e-13
    T=.018
    B=lambda x:1/np.expm1((x+.02)/T)
    F=lambda x:expit(-(x-.01)/T)
    bal=[]
    for p in (0.,3e-5,.012):
        coll,ev=collision(p,B,[],F,[],.82,(13,8,8,8),details=True)
        bal.append([ev,float(abs(coll[0]/coll[1]-1))])
    assert max(max(x) for x in bal)<4e-14
    return dict(independent_original_gamma_trace=trace_control(),differential_cross_transfer_defect=max(tests),
                thermal_event_and_integrated_balance_defects=bal)


def run(D):
    matches=read(D/"warm-near-rest-light-inverse-decay-2026-11-26.json.gz")
    res=[]
    for q in (16,32):
        saved=read(D/f"warm-light-occupation-next-q{q}-tape-2026-11-26.json.gz")
        saved+=[j["checkpoint"] for j in read(D/f"warm-observed-suboperator-q{q}-tape-2026-11-26.json.gz")["blocks"]]
        parent=next(a["stages"] for a in matches["cases"] if a["q"]==q)
        rows=[]
        for s,d in zip(saved,parent,strict=True):
            S=s["snapshot"]; h=s["snapshot_sha256"]
            assert h==digest(S)==d["snapshot_sha256"]
            B,be,F,fe=snapshot_functions(S)
            H=d["fine"]["H_over_rest_decay_width"]*width
            rest=[]
            for eta in (.7,.82,.9):
                rates=[]
                for n in (8,14,22):
                    a,_=collision(0.,B,be,F,fe,eta,(n,1,n,1))
                    rates.append((a/H).tolist())
                rate=np.asarray(rates)
                net=rate[:,0]-rate[:,1]
                assert np.all(net<0)
                assert abs(net[-1]/net[-2]-1)<.012
                rest.append(dict(sqrt_s_over_M_cut=eta,orders=[8,14,22],gain_loss_per_Hubble=rates,
                                 net_per_Hubble=net.tolist()))
            rows.append(dict(T_GeV=d["T_GeV"],snapshot_sha256=h,rest=rest,
                             archived_decay_net_at_rest_per_Hubble=(lambda k:k[0]+k[1]-k[2])(
                                 d["fine"]["spontaneous_stimulated_and_inverse_at_rest_per_Hubble"])))
        res.append(dict(q=q,stages=rows))
    return dict(scope="Independent frozen exact-mass light plus scalar elastic rates below fixed "
                      "CM energy cut. No pole/RIS or above-cut scattering, other elastic sources, "
                      "or full dynamical sign inferred.",controls=controls(),cases=res)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tapes",type=Path,default=Path("results"))
    print(json.dumps(run(p.parse_args().tapes),indent=2,allow_nan=False))
