"""Uniform scalar-pair sign on a declared time interpolation of accepted shapes."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import expit

from three_species_decay_onsager_reference import M,m,mp,y
from warm_light_occupation_stage_audit_reference import digest,read


def panels(p,P,a,b,mass,n):
    """Interval bounds as p/a decreases and entropy interpolates bilinearly."""
    # Earlier p, later P. u is the within-shell fraction in kinetic energy.
    u=np.linspace(0.,1.,n+1)
    def edges(x):
        e=np.hypot(x,mass)
        return e[:-1,None]+np.diff(e)[:,None]*u
    early,late=edges(p),edges(P)
    a,b=np.array(a),np.array(b)
    z1=a[0,:,None]+a[1,:,None]*u
    z2=b[0,:,None]+b[1,:,None]*u
    zz=np.stack((z1[:,:-1],z1[:,1:],z2[:,:-1],z2[:,1:]))
    return late[:,:-1].ravel(),early[:,1:].ravel(),zz.min(axis=0).ravel(),zz.max(axis=0).ravel()


def bounds(S,T,n=4):
    p,P=[np.asarray(s["lab_momentum_edges_GeV"]) for s in (S,T)]
    np.testing.assert_allclose(p[1:]/P[1:],
                               S["scale_factor"]**-1*T["scale_factor"],rtol=1.5e-13)
    assert np.all(P<=p) and T["scale_factor"]>S["scale_factor"]
    shapes=lambda x,b: x["entropy_variable_"+b]
    B=panels(p,P,shapes(S,"scalar"),shapes(T,"scalar"),mp,n)
    F=panels(p,P,shapes(S,"fermion"),shapes(T,"fermion"),m,n)
    Emin,Emax,zmin,zmax=B
    fl,fu,_,fz=F
    L=Emin[:,None]+Emin
    U=Emax[:,None]+Emax
    Z=zmin[:,None]+zmin
    et_hi=np.hypot(m,p[-1])
    et_lo=np.hypot(m,P[-1])
    eps=2e-14
    gap=np.inf
    for l,u,z in zip(fl,fu,fz,strict=True):
        ok=(L-u<=et_hi+eps)&(U-l>=m-eps)
        if np.any(ok):
            gap=min(gap,np.min(Z[ok])-z)
    assert np.isfinite(gap)
    if gap<=0.:
        return dict(subdivisions=n,represented_gap=float(gap),successful=False)

    F0,F1=np.hypot(m,p),np.hypot(m,P)
    s0=2*m*m+2*m*F1[:-1]
    s1=2*m*m+2*(et_hi*F0[1:]+p[-1]*p[1:])
    a=M*M-m*m-mp*mp+s1/2
    amp=8*y**4*np.maximum(0.,s0/4-m*m)*M*M/a**2
    flux=.5*np.sqrt(1-4*mp*mp/s0)
    fm=expit(-fz.reshape(len(p)-1,n).max(axis=1))
    negative=float(np.sum(amp*flux*fm*P[:-1]*np.diff(F1))/(128*np.pi**3*et_hi)*(-np.expm1(-gap)))
    Bmax=np.hypot(mp,p[-1])
    emax=2*Bmax-m
    if emax<=et_lo:
        upper=0.
    else:
        allow=Emax>=m+et_lo-Bmax-eps
        bmin=np.min(zmin[allow])
        select=(U>=m+et_lo-eps)&allow[:,None]&allow
        bsum=np.min(Z[select])
        bosons=np.exp(-bsum)/(-np.expm1(-bmin))**2
        smax=(2*Bmax)**2
        pp=np.sqrt(smax/4-m*m)
        A=M*M-m*m-mp*mp
        amp=2*y**4*((smax/4)*(2*pp/A)**2+(pp*(4*m+2*M)/A)**2)
        upper=float(amp*bosons*(emax-et_lo)*np.sqrt((emax-m)*(emax+m))/(128*np.pi**3*m))
    return dict(subdivisions=n,represented_gap=float(gap),negative_lower_GeV=negative,
                positive_upper_GeV=upper,relative_upper=upper/negative,successful=upper<negative)


def run(seed,work,k):
    base=read(seed)
    r=base[-1] if isinstance(base,list) else base
    out=[]
    for i in range(1,k+1):
        q=read(work/f"step-{i}.json")
        assert q.get("previous",q.get("parent_snapshot_sha256"))==r["snapshot_sha256"]
        assert digest(r["snapshot"])==r["snapshot_sha256"] and digest(q["snapshot"])==q["snapshot_sha256"]
        a=bounds(r["snapshot"],q["snapshot"],4)
        if not a["successful"]:
            a=bounds(r["snapshot"],q["snapshot"],8)
        out.append(dict(step=i,T_GeV=[r["snapshot"]["T_visible_GeV"],q["snapshot"]["T_visible_GeV"]],
                        snapshot_sha256=[r["snapshot_sha256"],q["snapshot_sha256"]],**a))
        r=q
    return dict(scope="Only the prescribed bilinear-in-time and shell-kinetic-coordinate entropy interpolation is bounded. Actual intermediate solutions need not follow it.",
                accepted_brackets=k,successful_brackets=sum(x["successful"] for x in out),brackets=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("seed",type=Path)
    p.add_argument("work",type=Path)
    p.add_argument("steps",type=int)
    a=p.parse_args()
    print(json.dumps(run(a.seed,a.work,a.steps),indent=2,allow_nan=False))
