"""Conditional finite joint comparison of absolute cold response profiles."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

from early_pole_endpoint_cold_normalized_profile_envelope import fit_state, geometry
from early_pole_endpoint_cold_common_moment_profile_test import vol
from early_pole_endpoint_cold_panel_reweight_envelope import mean_vertices, panel_excess


FIELDS = ((0, 0), (0, 1), (2, 0))
NAMES = ('forward_number', 'inverse_number', 'forward_excess')


def make(e, z, state, bins, shared):
    d, f, mass, _ = state
    a, b = e[:-1], e[1:]
    mid = (a+b)/2
    i = np.searchsorted(d.ds.ed, mid)-1
    j = np.searchsorted(bins, mid)-1
    if (np.any(i < 0) or np.any(i >= d.ns) or np.any(j < 0) or np.any(j >= len(bins)-1)
            or np.any(a < d.ds.ed[i]) or np.any(b > d.ds.ed[i+1])
            or np.any(a < bins[j]) or np.any(b > bins[j+1])):
        raise ArithmeticError('Unaligned panel')
    V = vol(a, b)
    O = np.zeros((len(FIELDS), len(bins)-1))
    H = np.zeros((len(a), *O.shape))
    original = np.zeros_like(O)
    for k, (race, sign) in enumerate(FIELDS):
        values = z[:, race, sign]
        s = f[i]
        if sign and np.any((s == 0) & (values != 0)):
            raise ArithmeticError('Inverse reaction in a zero occupation cell')
        if sign and np.any(s == 0) and shared is not None:
            ce, cf, _, *_ = shared
            idx = np.searchsorted(ce, mid)-1
            if np.any((s == 0) & (cf[idx]>0)):
                raise ArithmeticError('No known inverse coefficient for a newly occupied panel')
        coef = (values/(1+s) if not sign else
                np.divide(values, s, where=s>0, out=np.zeros_like(s)))
        if not sign:
            np.add.at(O[k], j, coef)
        H[np.arange(len(a)), k, j] = coef/V
        np.add.at(original[k], j, values)
    if np.max(abs(O+np.tensordot(V*f[i], H, axes=1)-original)/np.maximum(original, 1e-90)) > 4e-14:
        raise ArithmeticError('Original rates do not reconstruct')
    if shared is None:
        ce = np.r_[a[0], d.ds.ed[(d.ds.ed>a[0]) & (d.ds.ed<b[-1])], b[-1]]
        cf = f[np.searchsorted(d.ds.ed, (ce[:-1]+ce[1:])/2)-1]
    else:
        ce, cf, mass, *fixed_means = shared
    parent = np.searchsorted(ce, mid)-1
    if (np.any(parent<0) or np.any(parent>=len(cf)) or
            np.any(a<ce[parent]) or np.any(b>ce[parent+1])):
        raise ArithmeticError('No matching prescribed moment parent')
    X = panel_excess(a, b, mass, 48)
    pre = panel_excess(ce[:-1], ce[1:], mass, 48)
    means = fixed_means[0] if shared is not None and fixed_means else pre
    count, both = [], []
    for g, s in enumerate(cf):
        sel = parent == g
        target = vol(ce[g], ce[g+1])
        if abs(V[sel].sum()/target-1)>8e-14 or abs((V[sel] @ X[sel])/target/pre[g]-1)>8e-14:
            raise ArithmeticError('Conditioning count or mean quadrature failed')
        if s<=0:
            continue
        P = H[sel].reshape(np.count_nonzero(sel), -1)*target*s
        count.append(P)
        if len(P)==1 and abs(X[sel][0]/means[g]-1)>8e-14:
            raise ArithmeticError('Single panel does not attain prescribed mean')
        W = mean_vertices(X[sel], means[g]) if len(P)>1 else np.ones((1,1))
        both.append(W@P)
    return O, original, count, both


def gradients(card, member, ix, factors):
    start, _, *vertices = card
    out = np.concatenate(vertices[member], axis=0)[:, ix]/factors
    sizes = [len(z) for z in vertices[member]]
    E = np.zeros((len(sizes), sum(sizes)))
    for r, (p, q) in enumerate(zip(np.cumsum([0]+sizes)[:-1], np.cumsum(sizes), strict=True)):
        E[r, p:q] = 1.
    return start.ravel()[ix]/factors, out.T, E


def solve(a,b,field,member,*,method="highs",independent_bracket=False):
    """One common occupancy per chart must serve all requested fields."""
    totals=a[1].sum(axis=1)
    n=a[0].shape[1]
    if min(totals)<=0:
        raise ArithmeticError('Original left-chart reference totals must be positive')
    idx=np.concatenate([np.arange(i*n,(i+1)*n) for i in field])
    factor=np.repeat(totals[list(field)],n)
    U,F,P=gradients(a,member,idx,factor)
    V,G,Q=gradients(b,member,idx,factor)
    D=U-V
    k,l=F.shape[1],G.shape[1]
    j=len(D)
    delta=np.c_[F,-G]
    A=np.block([[delta,-np.eye(j)],[-delta,-np.eye(j)]])
    ub=np.r_[-D,D]
    E=np.block([[P,np.zeros((len(P),l+j))],
                [np.zeros((len(Q),k)),Q,np.zeros((len(Q),j))]])
    rhs=np.ones(len(E))
    c=np.r_[np.zeros(k+l),np.ones(j)]
    z=linprog(c,A_ub=A,b_ub=ub,A_eq=E,b_eq=rhs,bounds=(0.,None),
              method=method,options={'primal_feasibility_tolerance':1e-10,
              'dual_feasibility_tolerance':1e-10})
    if not z.success:
        raise ArithmeticError(z.message)
    x=z.x
    # Independently project each group of vertex weights to its simplex.
    # Report the direct objective at this concrete near-feasible witness.
    def projected(raw,A):
        out=np.maximum(raw,0.)
        scale=A@out
        if min(scale)<=0:
            raise ArithmeticError('Degenerate approximate primal witness')
        out/=scale@A
        return out
    xleft=projected(x[:k],P)
    xright=projected(x[k:k+l],Q)
    left=U+F@xleft
    right=V+G@xright
    pp,qq=z.eqlin.marginals,z.ineqlin.marginals
    # Audit a lower bound directly on the stored convex vertices. For any
    # |w| <= 1, w dot (left-right) <= ||left-right||_1. Its minimum separates
    # over the scalar parents. Clipping makes this independent of the solver's
    # feasibility errors in the dual multipliers.
    w=np.clip(qq[j:]-qq[:j],-1.,1.)
    ll=w@F; rr=-w@G
    m,n0=[np.array([min(t[row!=0]) for row in C]) for t,C in ((ll,P),(rr,Q))]
    support=float(np.longdouble(w)@np.longdouble(D)
                  +m.sum(dtype=np.longdouble)+n0.sum(dtype=np.longdouble))
    verb=(abs(left-right).reshape(len(field),n)).sum(axis=1)
    dual=rhs@pp+ub@qq
    residuals=[abs(z.fun-verb.sum()),float(np.max(abs(E@x-rhs))),float(max(0.,max(A@x-ub),-min(x))),
               float(max(0.,max(E.T@pp+A.T@qq-c),max(qq))),abs(z.fun-dual),abs(z.fun-support)]
    if independent_bracket:
        # The LP proposes two candidates, but accuracy is measured only by
        # actual projected simplex weights and a recomputed signed support.
        proj=max(float(np.max(abs(P@xleft-1))),float(np.max(abs(Q@xright-1))))
        hi=float(fsum(verb))
        lo=max(0.,support)
        if min(xleft)<0 or min(xright)<0 or proj>2e-14 or not -2e-12<=hi-lo<=2e-9:
            raise ArithmeticError('Independent finite lower/upper bracket did not resolve')
        return dict(minimum_l1=(hi+lo)/2,support_hyperplane_value=lo,projected_witness_value=hi,
                    verified_bracket_width=max(0.,hi-lo),projected_simplex_defect=proj,
                    optimizer_nominal_residuals=[float(x) for x in residuals],
                    support_signed_weights_by_field=w.reshape(len(field),n).tolist(),
                    witness_l1_by_field=verb.tolist())
    if max(residuals)>2e-9 or support>verb.sum()+2e-12:
        raise ArithmeticError(('Unverified optimization',residuals))
    return dict(minimum_l1=float(z.fun), witness_l1_by_field=verb.tolist(),
                field_signed_sum_delta_left_minus_right=(left-right).reshape(len(field),n).sum(axis=1).tolist(),
                verified_dual_value=float(dual), support_hyperplane_value=support,
                projected_witness_value=float(verb.sum()),
                support_signed_weights_by_field=w.reshape(len(field),n).tolist(),
                witness_abs_residuals_by_field=abs(left-right).reshape(len(field),n).tolist(),
                max_verification_defect=max(residuals))


def run(left,right,states,archive):
    charts=[fit_state(p,archive,n,raw) for p,n,raw in zip(states,(2,4),(left,right),strict=True)]
    d,y,m,_=charts[0]
    lo,hi=1e-6,1.5e-5
    ce=np.r_[lo,d.ds.ed[(d.ds.ed>lo)&(d.ds.ed<hi)],hi]
    share=ce,y[np.searchsorted(d.ds.ed,(ce[:-1]+ce[1:])/2)-1],m
    output=[]
    for t in (0,1):
        Z=list(map(geometry,(left[t],right[t])))
        bins=np.intersect1d(Z[0][0],Z[1][0])
        if len(bins)<5 or [bins[0],bins[-1]]!=[lo,hi]:
            raise ArithmeticError('Incomplete common support')
        data={}
        for label, shared in [('own',None),('shared_S2',share)]:
            a,b=[make(*z,s,bins,shared) for z,s in zip(Z,charts,strict=True)]
            for j,member in enumerate(('count_only','count_and_mean')):
                r={name:solve(a,b,(i,),j) for i,name in enumerate(NAMES)}
                r['joint']=solve(a,b,(0,1,2),j)
                x=np.array([r[name]['minimum_l1'] for name in NAMES])
                v=np.array(r['joint']['witness_l1_by_field'])
                if np.any(x>v+2e-9):
                    raise ArithmeticError('Joint residual violates its individual minima')
                data[label+'_'+member]=r
        orig=(np.abs(a[1]-b[1]).sum(axis=1)/a[1].sum(axis=1)).tolist()
        output.append((bins.tolist(),data,orig))
    (l,x,_),(r,y,orig)=output
    err=max(abs(x[k][j]['minimum_l1']-y[k][j]['minimum_l1']) for k in x for j in x[k])
    if l!=r or err>2e-9:
        raise ArithmeticError(('Meshes did not agree',err))
    return dict(x=.05,field_order=list(NAMES),common_q_edges=r,original_l1_by_field=orig,
                paired_mesh_optimum_maximum_defect=err,cases=y,
                scope='Conditional finite q panels and held chart-specific optical fits. Each field divided by its own original S2 band total, not each chart total. No directed rounding or physical bound.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('left','right','states'):
        p.add_argument('--'+k,nargs=2,type=Path,required=True)
    p.add_argument('--archive',type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.left,a.right,a.states,a.archive),indent=2,allow_nan=False))
