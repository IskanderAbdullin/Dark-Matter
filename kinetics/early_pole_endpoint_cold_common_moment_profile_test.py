"""Counterfactual finite S2/S4 profiles sharing the same coarser scalar moments."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_cold_normalized_profile_envelope import fit_state,geometry
from early_pole_endpoint_cold_panel_reweight_envelope import mean_vertices,panel_excess
from early_pole_endpoint_cold_profile_exact_finite_distance import solve


def vol(a,b):
    return (b-a)*(b*b+b*a+a*a)/3


def build(common,e0,reference,book,pop):
    r0,r1,counts,means,m=reference
    (e,z),(d,y,_,_)=book,pop
    a,b=e[:-1],e[1:]
    mid=(a+b)/2
    index=np.searchsorted(d.ds.ed,mid)-1
    parent=np.searchsorted(e0,mid)-1
    bins=np.searchsorted(common,mid)-1
    if (np.any(index<0) or np.any(index>=d.ns) or np.any(parent<0) or np.any(parent>=len(counts))
        or np.any(bins<0) or np.any(bins>=len(common)-1)
        or np.any(a<d.ds.ed[index]) or np.any(b>d.ds.ed[index+1])
        or np.any(a<r0[parent]*(1-5e-15)) or np.any(b>r1[parent]*(1+5e-15))
        or np.any(a<common[bins]) or np.any(b>common[bins+1])):
        raise ArithmeticError('A raw q panel straddles a common scalar parent or comparison bin')
    V=vol(a,b)
    A=z[:,0,0]/(1+y[index])
    C=np.bincount(bins,weights=A,minlength=len(common)-1)
    raw=np.bincount(bins,weights=z[:,0,0],minlength=len(common)-1)
    X=panel_excess(a,b,m,48)
    q=np.zeros((len(a),len(common)-1))
    q[np.arange(len(a)),bins]=A/V
    u,w=[],[]
    for j,(N,M) in enumerate(zip(counts,means,strict=True)):
        k=parent==j
        if abs(fsum(V[k])/vol(r0[j],r1[j])-1)>5e-14:
            raise ArithmeticError('A declared common scalar cell has not been covered')
        if not N>0:
            continue
        if abs(fsum(V[k]*X[k])/fsum(V[k])/M-1)>5e-14:
            raise ArithmeticError('Reference common mean is not feasible in both panel families')
        t=mean_vertices(X[k],M) if np.count_nonzero(k)>1 else np.ones((1,1))
        u.append(N*q[k])
        w.append(N*(t@q[k]))
    return dict(offset=C,total=fsum(raw),base=raw/fsum(raw),count_only=u,count_and_mean_excess=w)


def run(a,b,states,archive):
    charts=[fit_state(p,archive,n,raw) for p,n,raw in zip(states,(2,4),(a,b),strict=True)]
    old,y,m,_=charts[0]
    e0=np.union1d(old.ds.ed,[1e-6,1.5e-5])
    e0=e0[(e0>=1e-6)&(e0<=1.5e-5)]
    l,r=e0[:-1],e0[1:]
    j=np.searchsorted(old.ds.ed,(l+r)/2)-1
    reference=l,r,vol(l,r)*y[j],panel_excess(l,r,m,48),m
    out=[]
    for t in (0,1):
        P,Q=map(geometry,(a[t],b[t]))
        cut=np.intersect1d(P[0],Q[0])
        L,R=[build(cut,e0,reference,z,s) for z,s in zip((P,Q),charts,strict=True)]
        val=dict(common_q_edges=cut.tolist(),**{k:solve(L,R,k) for k in ('count_only','count_and_mean_excess')})
        out.append(val)
    lo,hi=out
    key='minimum_finite_common_bin_TV'
    err=max(abs(lo[m][key]-hi[m][key]) for m in ('count_only','count_and_mean_excess'))
    if err>2e-9:
        raise ArithmeticError('Common moment control disagrees between accepted meshes')
    return dict(x=.05,charge_div=2,reference_scalar_shells=2,fixed_excess_definition_mass=m,**hi,max_pair_mesh_difference=err,
                scope='Counterfactual common coarser moments: releases original S4 fine-shell constraints and allows each chart its own held optical fits. Existence within finite normalized profiles only, not actual evolved agreement.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('left','right','states'):
        p.add_argument('--'+k,nargs=2,type=Path,required=True)
    p.add_argument('--archive',type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.left,a.right,a.states,a.archive),indent=2,allow_nan=False))
