"""Test robustness to the reference used for common projected scalar moments."""
import argparse
import json
from pathlib import Path

import numpy as np

from early_pole_endpoint_cold_joint_absolute_profile import FIELDS, make, solve
from early_pole_endpoint_cold_normalized_profile_envelope import fit_state, geometry
from early_pole_endpoint_cold_common_moment_profile_test import vol
from early_pole_endpoint_cold_panel_reweight_envelope import panel_excess


def project(state, reference, m):
    d, f, _, _ = state
    lo, hi = reference[0], reference[-1]
    ce=np.union1d(reference, d.ds.ed)
    ce=ce[(lo<=ce)&(ce<=hi)]
    a,b=ce[:-1],ce[1:]
    mid=(a+b)/2
    i=np.searchsorted(d.ds.ed,mid)-1
    j=np.searchsorted(reference,mid)-1
    if (np.any(i<0) or np.any(i>=d.ns) or np.any(j<0) or np.any(j>=len(reference)-1)
            or np.any(a<d.ds.ed[i]) or np.any(b>d.ds.ed[i+1])
            or np.any(a<reference[j]) or np.any(b>reference[j+1])):
        raise ArithmeticError('Invalid exact overlap')
    V=vol(a,b)
    M=panel_excess(a,b,m,48)
    weighted=V*f[i]
    n=np.bincount(j,weights=weighted,minlength=len(reference)-1)
    u=np.bincount(j,weights=weighted*M,minlength=len(reference)-1)
    volume=np.bincount(j,weights=V,minlength=len(reference)-1)
    if np.max(abs(volume/vol(reference[:-1],reference[1:])-1))>3e-14:
        raise ArithmeticError('Moment projection loses phase-space volume')
    return n,u


def feasibility(e,ce,m,means,n):
    a,b=e[:-1],e[1:]
    idx=np.searchsorted(ce,(a+b)/2)-1
    if np.any(idx<0) or np.any(idx>=len(n)) or np.any(a<ce[idx]) or np.any(b>ce[idx+1]):
        raise ArithmeticError('A panel straddles the prescribed parents')
    x=panel_excess(a,b,m,48)
    errors=[]
    for i in np.flatnonzero(n>0):
        y=x[idx==i]
        if not len(y):
            raise ArithmeticError('Missing positive reference parent')
        mn,mx=min(y),max(y)
        if not mn*(1-8e-14)<=means[i]<=mx*(1+8e-14):
            errors.append(dict(parent_edges=ce[i:i+2].tolist(), requested_mean=means[i],
                               attainable_interval=[mn,mx],
                               relative_excess_distance=max(mn-means[i],means[i]-mx)/means[i]))
    return errors


def max_feasible_common_partition(ce, fields, charts, m):
    """Retain as many original conditioning boundaries as possible."""
    x=[]
    for e in charts:
        a,b=e[:-1],e[1:]
        parent=np.searchsorted(ce,(a+b)/2)-1
        if np.any(parent<0) or np.any(parent>=len(ce)-1) or np.any(a<ce[parent]) or np.any(b>ce[parent+1]):
            raise ArithmeticError('Incompatible parent candidate')
        x.append((parent,panel_excess(a,b,m,48)))
    best=[None]*len(ce)
    best[0]=[0]
    for j in range(1,len(ce)):
        for i in range(j):
            if best[i] is None:
                continue
            caps=[]
            for k,z in x:
                w=z[(i<=k)&(k<j)]
                if not len(w):
                    raise ArithmeticError('Missing complete parent region')
                caps.append((min(w),max(w)))
            low,high=max(t[0] for t in caps),min(t[1] for t in caps)
            means=[sum(u[i:j])/sum(n[i:j]) for n,u in fields if sum(n[i:j])>0]
            if all(low*(1-8e-14)<=y<=high*(1+8e-14) for y in means):
                trail=best[i]+[j]
                if best[j] is None or len(trail)>len(best[j]):
                    best[j]=trail
    if best[-1] is None:
        raise ArithmeticError('Even the whole-band projected moments have no common support')
    return np.asarray(best[-1])


def run(left,right,states,archive):
    charts=[fit_state(p,archive,n,raw) for p,n,raw in zip(states,(2,4),(left,right),strict=True)]
    d,f,m,_=charts[0]
    lo,hi=1e-6,1.5e-5
    ce=np.r_[lo,d.ds.ed[(d.ds.ed>lo)&(d.ds.ed<hi)],hi]
    V=vol(ce[:-1],ce[1:])
    (n0,u0),(n1,u1)=[project(s,ce,m) for s in charts]
    original=f[np.searchsorted(d.ds.ed,(ce[:-1]+ce[1:])/2)-1]
    if np.any(abs(n0-V*original)>5e-14*n0):
        raise ArithmeticError('Projected left moments disagree with stored left chart')
    alphas=np.linspace(0.,1.,5)
    meshes=[]
    for t in (0,1):
        Z=list(map(geometry,(left[t],right[t])))
        bins=np.intersect1d(Z[0][0],Z[1][0])
        if len(bins)<5 or [bins[0],bins[-1]]!=[lo,hi]:
            raise ArithmeticError('Incomplete common support')
        data=[]
        for alpha in alphas:
            N=(1-alpha)*n0+alpha*n1
            U=(1-alpha)*u0+alpha*u1
            means=np.divide(U,N,out=panel_excess(ce[:-1],ce[1:],m,48),where=N>0)
            fails=[feasibility(z[0],ce,m,means,N) for z in Z]
            row=dict(S4_reference_weight=float(alpha), mean_feasibility_failures_by_chart=fails)
            counts=[make(*z,s,bins,(ce,N/V,m)) for z,s in zip(Z,charts,strict=True)]
            row['shared_count_only']=solve(*counts,tuple(range(len(FIELDS))),0)
            if not any(fails):
                both=[make(*z,s,bins,(ce,N/V,m,means)) for z,s in zip(Z,charts,strict=True)]
                row['shared_count_and_mean']=solve(*both,tuple(range(len(FIELDS))),1)
            data.append(row)
        idx=max_feasible_common_partition(ce,[(n0,u0),(n1,u1)],[z[0] for z in Z],m)
        ne=ce[idx]
        NV=vol(ne[:-1],ne[1:])
        folded=[]
        for alpha in alphas:
            N=(1-alpha)*n0+alpha*n1
            U=(1-alpha)*u0+alpha*u1
            n=np.add.reduceat(N,idx[:-1])
            u=np.add.reduceat(U,idx[:-1])
            means=np.divide(u,n,out=panel_excess(ne[:-1],ne[1:],m,48),where=n>0)
            if any(feasibility(z[0],ne,m,means,n) for z in Z):
                raise ArithmeticError('Coarsened witness is not feasible')
            cards=[make(*z,s,bins,(ne,n/NV,m,means)) for z,s in zip(Z,charts,strict=True)]
            folded.append(dict(S4_reference_weight=float(alpha), both=solve(*cards,tuple(range(len(FIELDS))),1)))
        meshes.append((data,idx,folded))
    (a,ia,ca),(b,ib,cb)=meshes
    if not np.array_equal(ia,ib):
        raise ArithmeticError('Paired meshes yield different maximal refinements')
    defect=0.
    for p,q in zip(a,b,strict=True):
        if p.keys()!=q.keys() or p['mean_feasibility_failures_by_chart']!=q['mean_feasibility_failures_by_chart']:
            raise ArithmeticError('The paired meshes disagree about prescribed moment attainability')
        for field in p.keys()-{'S4_reference_weight','mean_feasibility_failures_by_chart'}:
            defect=max(defect,abs(p[field]['minimum_l1']-q[field]['minimum_l1']))
    defect=max(defect,max(abs(x['both']['minimum_l1']-y['both']['minimum_l1']) for x,y in zip(ca,cb,strict=True)))
    if defect>2e-9:
        raise ArithmeticError('Inconsistent absolute response comparison')
    return dict(x=.05, fixed_reference_mass=m,coarse_parent_edges=ce.tolist(),
                reference_count=[n0.tolist(),n1.tolist()],reference_excess=[u0.tolist(),u1.tolist()],
                maximum_optimum_difference=defect, rows=b,
                common_maximal_partition_indices=ib.tolist(), coarsened_rows=cb,
                scope='Finite convex mixtures of the two projected moment vectors, not fits or perturbations of either physical evolution.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('left','right','states'):
        p.add_argument('--'+k,nargs=2,type=Path,required=True)
    p.add_argument('--archive',type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.left,a.right,a.states,a.archive),indent=2,allow_nan=False))
