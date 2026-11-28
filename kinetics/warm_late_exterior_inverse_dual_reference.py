"""Use a per-cohort positive affine majorant for mixed inverse reactions."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import expit, roots_legendre

from joint_visible_heavy_adiabatic_reference import snapshot_functions

from factorized_prompt_return_graph_reference import C
from frw_exterior_restart_reference import import_restart
from late_radiation_background_reference import MPL,TabulatedRadiation
from three_species_decay_onsager_reference import M,m,mp
from warm_late_mixed_inverse_moment_envelope_reference import evaluate
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom


def panel(S,cell,cut,group,n):
    ds,de=(mp,m) if cell==0 else (m,mp)
    g=4 if cell==0 else 1
    p=np.asarray(S["lab_momentum_edges_GeV"])
    ep=np.hypot(ds,p)
    key="entropy_variable_"+("scalar" if cell==0 else "fermion")
    a,b=np.asarray(S[key])
    s=a+np.minimum(b,0.)
    sup=np.exp(-s)/(-np.expm1(-s)) if cell==0 else expit(-s)
    x,d,_,side,floor,N,mean=group
    assert np.all((d==de)&(side==1))
    L=np.minimum(floor,mean)
    z=de*np.sqrt((1-x)*(1+x))
    current_floor=np.hypot(z,x*L)
    lo0=np.min(current_floor)
    margins=np.clip(cut-ep,lo0,8.)
    grid=np.unique(np.r_[np.geomspace(lo0,8.,n+1),margins,current_floor,np.inf])
    lo,up=grid[:-1],grid[1:]
    upper=np.isfinite(up)
    D=M*M-ds*ds-de*de
    r=np.sqrt((D-2*ds*de)*(D+2*ds*de))
    atmin=np.clip(D/(2*ds),lo,up)
    q=np.sqrt((atmin-de)*(atmin+de))
    species_low=(D*D+4*ds*ds*q*q)/(2*(D*atmin+r*q))
    v=up[upper]
    species_high=np.full_like(up,np.inf)
    species_high[upper]=(D*v+r*np.sqrt((v-de)*(v+de)))/(2*de*de)
    low=np.maximum(ep[:-1],np.maximum(cut-up,species_low*(1.-2e-14))[:,None])
    high=np.minimum(ep[1:],(species_high*(1.+2e-14))[:,None])
    span=np.maximum(high-low,0.)
    I0=(span*sup).sum(axis=1)
    I1=(span*np.where(span>0.,(low+high)/2,0.)*sup).sum(axis=1)
    U=np.sqrt((lo-de)*(lo+de))
    # Bound (E+x)/(p_ext E) x-by-x. The final semi-infinite
    # panel uses the same bound, not an arbitrary high momentum cap.
    Q=(I0+I1/lo)/U
    wedge=np.zeros_like(I0)
    wedge[upper]=np.maximum(0.,I1[upper]+(up[upper]-cut)*I0[upper])
    wedge[~upper]=I1[~upper]
    W=np.minimum(np.minimum(I1,cut*I0),wedge)/(lo*U)
    if n==48:
        v=np.geomspace(lo0,1e5,143)
        i=np.minimum(np.searchsorted(grid,v,side="right")-1,len(lo)-1)
        p1=np.sqrt((v-de)*(v+de))
        cmin=(D*D+4*ds*ds*p1*p1)/(2*(D*v+r*p1))
        cmax=(D*v+r*p1)/(2*de*de)
        f,s0,gfun,s1=snapshot_functions(S)
        fn=(f,gfun)[cell]
        zl,zw=roots_legendre(16)
        left=np.maximum(ep[:-1],np.maximum(cut-v,cmin)[:,None])
        right=np.minimum(ep[1:],cmax[:,None])
        span=np.maximum(right-left,0.)
        sample=left[:,:,None]+span[:,:,None]*(1+zl)/2
        wgt=span[:,:,None]*zw/2*fn(sample)/(v*p1)[:,None,None]
        E=v[:,None,None]
        truth=(wgt*(E+sample)).sum(axis=(1,2))
        cost=np.minimum(np.minimum(sample,E),np.minimum(cut,np.maximum(E+sample-cut,0.)))
        growth=(wgt*cost).sum(axis=(1,2))
        assert np.all(truth<=Q[i]*(1+1e-12)+1e-240)
        assert np.all(growth<=W[i]*(1+1e-12)+1e-240)
    assert np.all(np.isfinite(Q)) and np.all(np.isfinite(W))
    eligible=up[None,:]>=current_floor[:,None]
    first=np.maximum(lo,current_floor[:,None])
    birth=np.sqrt((first-z[:,None])*(first+z[:,None]))/x[:,None]
    birth=np.maximum(0.,birth-L[:,None])
    spend=mean-L
    vals=np.stack((Q,W))
    out=[]
    for value in vals:
        V=np.broadcast_to(value,birth.shape)
        peak=np.max(np.where(eligible,V,0.),axis=1)
        base=peak/np.maximum(spend,1e-20)
        best=peak.copy()
        # Each beta produces a proved affine majorant. The finite
        # search tightens it but need not find the optimum.
        for f in np.geomspace(1e-5,1e2,85):
            beta=base*f
            alpha=np.max(np.where(eligible,V-beta[:,None]*birth,-np.inf),axis=1)
            best=np.minimum(best,alpha+beta*spend)
        assert np.min(best)>=0.
        out.append(float(C*2*np.pi**2/g*(N@best)))
    return out


def run(eos,prior,archive):
    q=read(archive)
    S=q["snapshot"]
    assert q["snapshot_sha256"]==digest(S)
    root=read(prior)["branches"][-1]
    cut=np.hypot(M,1./root["entropy_expansion_from_original_anchor"])
    ext=import_restart(S["positive_exterior_restart"])
    arr=ext.array()
    p=np.asarray(S["lab_momentum_edges_GeV"])
    rho=get_mom(p,[S["entropy_variable_"+j] for j in ("scalar","fermion")])[2].sum()
    H=np.sqrt((TabulatedRadiation(eos,S["T_visible_GeV"]).state(1.)[1]+rho+ext.energy(1.))/3)/MPL
    B=evaluate(S,eos,cut,np.inf,True)["controls"][-1]
    old=np.stack((B["normalized_bare_inverse_sector_matrix"],B["normalized_positive_excess_source_sector_matrix"]))
    rows=[]
    for n in (48,96,192,384):
        Mat=old.copy()
        for i,de in enumerate((m,mp)):
            new=np.asarray(panel(S,i,cut,arr[:,arr[1]==de],n))/(rho*H)
            Mat[:,i,1-i]=np.minimum(new,old[:,i,1-i])
        rows.append(dict(outer_positive_panels=n,bare_and_excess_matrix=Mat.tolist(),
                         selected_upper=Mat.sum(axis=(1,2)).tolist()))
    return dict(scope="Instantaneous upper for archived core occupations capped by their observed cell extrema "
                      "and unknown distributions with the preserved external cohort moments. No bound on unknown "
                      "prior sources, time evolution, or repeated reuse. The tiny exterior/exterior term uses "
                      "the independent group-moment envelope.",
                snapshot_sha256=q["snapshot_sha256"],threshold_energy_GeV=cut,
                reference_group_moment_upper=old.sum(axis=(1,2)).tolist(),controls=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","prior","archive"):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.prior,a.archive),indent=2,allow_nan=False))
