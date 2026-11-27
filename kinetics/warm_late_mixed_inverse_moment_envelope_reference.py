"""Bound the late high inverse rate with resolved moments and admitted cohorts."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import expit

from factorized_prompt_return_graph_reference import C
from frw_exterior_restart_reference import import_restart
from late_radiation_background_reference import MPL,TabulatedRadiation
from three_species_decay_onsager_reference import M,m,mp
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_high_parent_ledger_reference import high
from warm_prompt_reentry_shape_audit_reference import get_mom


def reach(a,b,la,ua,lb,ub,cut,use_angle):
    x,y=la.copy(),lb.copy()
    D=M*M-a*a-b*b
    gamma=np.sqrt((D-2*a*b)*(D+2*a*b))
    lb=np.maximum(lb,np.maximum(b,cut-ua))
    la=np.maximum(la,np.maximum(a,cut-ub))
    x=np.maximum(x,la)
    y=np.maximum(y,lb)

    def lower(d,other,z):
        q=np.sqrt(np.maximum((z-other)*(z+other),0.))
        return (D*D+4*d*d*q*q)/(2*(D*z+gamma*q))

    if use_angle:
        # Lower daughter E on the exact inverted angular ellipse as
        # a function of its partner E. It is convex and minimized at
        # partner E=D/(2*d); clip this first to necessary parent support.
        xv=np.maximum(lb,np.minimum(ub,D/(2*a)))
        yv=np.maximum(la,np.minimum(ua,D/(2*b)))
        x=np.maximum(x,lower(a,b,xv)*(1.-2e-14))
        y=np.maximum(y,lower(b,a,yv)*(1.-2e-14))
    return x,y,(x<=ua)&(y<=ub)&(ua+ub>=cut)


def sanity():
    rng=np.random.default_rng(2080)
    D=M*M-m*m-mp*mp
    passed=0
    for i in range(140):
        x=np.sort(np.r_[mp,rng.uniform(1e-4,1.8,60)])
        y=np.sort(np.r_[m,rng.uniform(.0201,1.8,65)])
        la,ua=x[:-1,None],x[1:,None]
        lb,ub=y[None,:-1],y[None,1:]
        a=la+(ua-la)*rng.random((len(x)-1,len(y)-1))
        b=lb+(ub-lb)*rng.random(a.shape)
        p=np.sqrt((a-mp)*(a+mp));q=np.sqrt((b-m)*(b+m))
        cut=rng.uniform(.4,1.)
        ok=(a+b>=cut)&(abs(a*b-D/2)<=p*q)
        u,v,s=reach(mp,m,la,ua,lb,ub,cut,True)
        assert np.all(s[ok]) and np.all(u[ok]<=a[ok]) and np.all(v[ok]<=b[ok])
        passed+=int(ok.sum())
    assert passed>100
    return passed


def species(S,cap):
    p=np.asarray(S["lab_momentum_edges_GeV"])
    theta=[np.asarray(S["entropy_variable_"+z]) for z in ("scalar","fermion")]
    N,K,E=get_mom(p,theta)
    ex=import_restart(S["positive_exterior_restart"]).array()
    out=[]
    original=len(p)-1
    for j,d in enumerate((mp,m)):
        ec=np.hypot(p,d)
        a,mass,_,side,bd,n,mean=ex[:,ex[1]==d]
        assert np.all(side==1)
        upper=np.hypot(cap,d)
        lower=np.minimum(bd,mean)
        lowerp=np.hypot(d*np.sqrt((1-a)*(1+a)),a*lower)
        assert np.all(np.hypot(d*np.sqrt((1-a)*(1+a)),a*mean)<=upper)
        ar=np.r_[np.ones(original),a]
        floor=np.r_[ec[:-1],lowerp]
        capv=np.r_[ec[1:],np.full(len(n),upper)]
        base=np.r_[ec[:-1],lower]
        excess=np.r_[np.maximum(0.,K[j]-N[j]*(ec[:-1]-d)),n*(mean-lower)]
        number=np.r_[N[j],n]
        a1,b1=theta[j]
        low=a1+np.minimum(0.,b1)
        if j==0:
            assert np.min(low)>0.
            ceiling=np.exp(-low)/(-np.expm1(-low))
        else:
            ceiling=expit(-low)
        out.append((ar,floor,capv,base,excess,number,np.r_[ceiling,np.full(len(n),np.inf)]))
    return out,original,float(E.sum())


def integral(mass,deg,group,t,eligible,use_moments,use_occupancy):
    ar,low,up,base,excess,N,sup=group
    n=np.broadcast_to(N,t.shape).copy()
    t=np.where(eligible,t,np.broadcast_to(low,t.shape))
    p=np.sqrt(np.maximum((t-mass)*(t+mass),0.))
    # Zero leg moment is innocuous. For a physically populated pair,
    # the high-parent and angular gates exclude the mass endpoint.
    assert not np.any(eligible&(p==0.)&(n>0.)&(~np.isfinite(sup) if use_occupancy else True))
    if use_moments:
        rest=mass*np.sqrt((1.-ar)*(1.+ar))
        required=np.sqrt(np.maximum((t-rest)*(t+rest),0.))/ar
        delta=required-base
        n=np.minimum(n,np.divide(excess,delta,out=n.copy(),where=delta>2e-14*np.maximum(base,1e-60)))
    scale=np.where(eligible&(p==0.)&(n>0.),np.inf,0.)
    I=np.divide(2*np.pi**2*n/deg,p*t,out=scale,where=eligible&(p>0.))
    J=I*t
    if use_occupancy:
        finite=np.broadcast_to(np.isfinite(sup),t.shape)&eligible
        cap=np.broadcast_to(up,t.shape)[finite]
        # This optional ceiling is on retained per-mode occupation only.
        # The unknown archived exterior is not assumed to have such a cap.
        edge=np.broadcast_to(sup,t.shape)[finite]
        v=np.maximum(cap-t[finite],0.)
        I[finite]=np.minimum(I[finite],edge*v)
        J[finite]=np.minimum(J[finite],edge*v*(cap+t[finite])/2)
    assert np.isfinite(I).all() and np.isfinite(J).all()
    return I,J


def evaluate(S,eos,cut,cap,use_occupancy=False):
    groups,k,rho=species(S,cap)
    a,b=[[v[:,None] if j==0 else v[None,:] for v in z] for j,z in enumerate(groups)]
    ext=import_restart(S["positive_exterior_restart"])
    H=np.sqrt((TabulatedRadiation(eos,S["T_visible_GeV"]).state(1.)[1]+rho+ext.energy(1.))/3)/MPL
    report=[]
    for use_m in (False,True):
        x,y,ok=reach(mp,m,a[1],a[2],b[1],b[2],cut,True)
        I,J=integral(mp,1,a,x,ok,use_m,use_occupancy)
        V,W=integral(m,4,b,y,ok,use_m,use_occupancy)
        raw=C*(I*W+J*V)/(rho*H)
        plus=C*np.minimum(np.minimum(J*V,I*W),cut*I*V)/(rho*H)
        ix=(slice(0,k),slice(k,None))
        R=np.array([[raw[i,j].sum() for j in ix] for i in ix])
        P=np.array([[plus[i,j].sum() for j in ix] for i in ix])
        report.append(dict(use_number_and_birth_or_cell_energy=use_m,
                           normalized_bare_inverse_sector_matrix=R.tolist(),
                           normalized_positive_excess_source_sector_matrix=P.tolist(),
                           normalized_bare_inverse_sum=float(R.sum()),
                           normalized_positive_excess_sum=float(P.sum())))
    assert report[1]["normalized_bare_inverse_sum"]<=report[0]["normalized_bare_inverse_sum"]
    native=high(S,np.sqrt((cut-M)*(cut+M)),10)[1]/(rho*H)
    assert 0.<native<=report[-1]["normalized_bare_inverse_sector_matrix"][0][0]
    return dict(admitted_momentum_cap_GeV=cap if np.isfinite(cap) else None,
                separate_retained_occupancy_ceiling=use_occupancy,
                native_retained_inverse_check=native,controls=report)


def run(eos,prior,archive):
    q=read(archive)
    S=q["snapshot"]
    assert digest(S)==q["snapshot_sha256"]
    x=read(prior)["branches"][-1]
    pcut=1./x["entropy_expansion_from_original_anchor"]
    data=[evaluate(S,eos,np.hypot(M,pcut),cap) for cap in (.70,.85,1.05,1.50,2.0)]
    unlimited=evaluate(S,eos,np.hypot(M,pcut),np.inf,True)
    return dict(scope="Instantaneous, shape-independent two-moment bound only for specified represented "
                      "shells and ballistic saved cohorts. Finite cases impose admitted support and no phase ceilings. "
                      "The alternative with unbounded exterior uses a separate resolved pointwise ceiling. "
                      "Does not bound rates of unarchived matter, later re-entry, or integrated evolution.",
                snapshot_sha256=q["snapshot_sha256"],parent_momentum_cut_GeV=pcut,
                core_or_exterior_indices=["retained","known exterior"],sampled_kinematically_valid_rectangles=sanity(),
                cases=data,unbounded_exterior_with_resolved_ceiling=unlimited)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for x in ("eos","prior","archive"):
        p.add_argument(x,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.prior,a.archive),indent=2,allow_nan=False))
