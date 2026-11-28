"""Lower bound q transport under both finite occupation families and arbitrary panel support."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.optimize import linprog

from early_pole_endpoint_cold_normalized_profile_envelope import build,fit_state,geometry
from early_pole_endpoint_cold_profile_exact_finite_distance import (
    as_system, normalized_support_lower, profile_from_projected_blocks,
)


def repair_flow(candidate,p,q):
    F=np.maximum(candidate,0.).copy()
    for G,target in ((F,p),(F.T,q)):
        for v,cap in zip(G,target,strict=True):
            mass=fsum(v)
            if mass>cap:
                v*=cap/mass
    a=np.maximum(0.,p-F.sum(axis=1))
    b=np.maximum(0.,q-F.sum(axis=0))
    N=fsum(a)
    K=fsum(b)
    if abs(N-K)>3e-14:
        raise ArithmeticError('Repaired transport residual masses do not agree')
    if N>0. and K>0.:
        F+=np.outer(a/N,b)
    err=max(float(np.max(abs(F.sum(axis=1)-p))),float(np.max(abs(F.sum(axis=0)-q))))
    if min(F.ravel())<0. or err>4e-14:
        raise ArithmeticError('Projected transport does not match both independent marginals')
    return F,err


def certify(left,right,mode,lp,rp,prior,sl,sr,cost):
    p,q=[profile_from_projected_blocks(c,mode,z) for c,z in ((left,lp),(right,rp))]
    correction=max(0.,float(np.max(sl[:,None]+sr[None,:]-cost)))+64*np.finfo(float).eps
    sl=sl-correction
    if np.any(sl[:,None]+sr[None,:]>cost):
        raise ArithmeticError('A repaired coupling potential exceeds local cost')
    guard=64*np.finfo(float).eps
    bot=max(0.,normalized_support_lower(left,mode,sl)+normalized_support_lower(right,mode,sr)-guard)
    flow,err=repair_flow(prior,p,q)
    top=float(np.sum(flow*cost))+guard
    if not bot<=top or top-bot>2e-9:
        raise ArithmeticError(('Independent interval transport bracket is too wide',bot,top))
    return bot,top,flow,err,p,q


def transport(left,right,ea,eb,mode,*,independent_bracket=False):
    p,P,u=as_system(left,mode)
    q,Q,v=as_system(right,mode)
    n,a=p.shape
    m,b=q.shape
    if len(ea)!=n+1 or len(eb)!=m+1 or np.min(np.diff(ea))<=0 or np.min(np.diff(eb))<=0:
        raise ArithmeticError('Require exact ordered support for every marginal')
    if ea[0]!=eb[0] or ea[-1]!=eb[-1]:
        raise ArithmeticError('Cold test support is not shared')
    unit=ea[-1]-ea[0]
    # Pay only the distance between the two intervals: adjacent or
    # overlapping panels may transport at zero, by design.
    cost=np.maximum(0.,np.maximum(eb[:-1][None,:]-ea[1:,None],
                                  ea[:-1,None]-eb[1:][None,:]))/unit
    mid=len(cost.ravel())
    j=np.arange(mid)
    row=sparse.coo_matrix((np.ones(mid),(np.repeat(np.arange(n),m),j)),shape=(n,mid))
    col=sparse.coo_matrix((np.ones(mid),(np.tile(np.arange(m),n),j)),shape=(m,mid))
    C=sparse.bmat([[sparse.csr_matrix(P),None,None],
                   [None,sparse.csr_matrix(Q),None],
                   [sparse.csr_matrix(-p),None,row],
                   [None,sparse.csr_matrix(-q[:-1]),col.tocsr()[:-1]]],format='csr')
    bb=np.r_[u,v,np.zeros(n+m-1)]
    z=np.r_[np.zeros(a+b),cost.ravel()]
    if independent_bracket:
        failures=[]
        for scale,method in ((256.,'highs'),(1024.,'highs-ipm'),(1.,'highs-ds')):
            G=linprog(z,A_eq=C,b_eq=bb*scale,bounds=(0.,None),method=method,
                      options={'primal_feasibility_tolerance':1e-10,'dual_feasibility_tolerance':1e-10})
            if not G.success:
                failures.append(G.message)
                continue
            w=G.x/scale
            support=G.eqlin.marginals
            ca=support[-(n+m-1):-(m-1)] if m>1 else support[-n:]
            cb=np.r_[support[-(m-1):],0.] if m>1 else np.zeros(1)
            try:
                L,H,T,err,pa,pb=certify(left,right,mode,w[:a],w[a:a+b],w[a+b:].reshape(n,m),ca,cb,cost)
            except ArithmeticError as e:
                failures.append(str(e))
                continue
            return dict(minimum_interval_transport_cost_q=(L+H)*unit/2.,
                        independent_cost_lower_upper_q=[L*unit,H*unit],
                        independent_cost_bracket_width_q=(H-L)*unit,
                        normalized_projected_marginal_max_absolute_defect=err,
                        left_right_profile_at_optimum=[pa.tolist(),pb.tolist()],
                        positive_transport_edges=int(np.count_nonzero(T)),
                        nominal_solver_constraint_defect=float(np.max(abs(C@w-bb))),
                        optimization_rhs_scale=scale)
        raise ArithmeticError(('No independently verified interval transport candidates',failures))
    g=linprog(z,A_eq=C,b_eq=bb,bounds=(0.,None),method='highs',
              options={'primal_feasibility_tolerance':1e-10,'dual_feasibility_tolerance':1e-10,'presolve':False})
    if not g.success:
        raise ArithmeticError(('Interval transport solve failed',g.message))
    w=g.x
    a0,b0=p@w[:a],q@w[a:a+b]
    flow=w[a+b:].reshape(n,m)
    err=max(np.max(abs(C@w-bb)),np.max(abs(flow.sum(axis=1)-a0)),np.max(abs(flow.sum(axis=0)-b0)),max(0.,-np.min(w)))
    lam=g.eqlin.marginals
    dual=max(float(np.max(C.T@lam-z)),0.)
    val=float(np.sum(flow*cost))
    lower=float(bb@lam)
    dg=max(abs(g.fun-lower),abs(g.fun-val))
    if max(err,dual,dg)>2e-9:
        raise ArithmeticError(('Interval transport primal and dual numerics did not close',err,dual,dg,g.fun))
    return dict(minimum_interval_transport_cost_q=float(g.fun*unit),
                dual_optimum_cost_q=lower*unit,realized_coupling_cost_q=val*unit,
                normalized_marginal_max_absolute_defect=float(err),
                normalized_dual_max_positive_defect=dual,normalized_cost_closure=dg,
                left_right_profile_at_optimum=[a0.tolist(),b0.tolist()],
                dimensionless_fractional_total_relative_to_original=[float(1/w[0]),float(1/w[a])],
                positive_transport_edges=int(np.count_nonzero(flow)))


def fixed(a):
    p=np.asarray(a['base'])
    return dict(offset=p,base=p,total=1.,fixed=[])


def run(raw,states,archive,levels,observable='forward_number',*,independent_bracket=False):
    if tuple(levels) not in ((2,4),(4,8),(8,16)):
        raise ArithmeticError('Use two accepted adjacent scalar levels')
    entries=[fit_state(p,archive,n,r) for p,n,r in zip(states,levels,raw,strict=True)]
    out=[]
    for j in (0,1):
        a,b=[geometry(p[j]) for p in raw]
        L,R=[build(g,t,g[0],raw=True,observable=observable) for g,t in zip((a,b),entries,strict=True)]
        e,f=a[0],b[0]
        base=transport(fixed(L),fixed(R),e,f,'fixed',independent_bracket=independent_bracket)
        U=transport(L,R,e,f,'count_only',independent_bracket=independent_bracket)
        V=transport(L,R,e,f,'count_and_mean_excess',independent_bracket=independent_bracket)
        k='minimum_interval_transport_cost_q'
        tol=5e-9*(e[-1]-e[0])
        if not U[k]-tol<=V[k]<=base[k]+tol:
            raise ArithmeticError('Nested constraints must not move the optimum in the wrong order')
        out.append(dict(unchanged_panel_amounts=base,count_only=U,count_and_mean_excess=V))
    err=max(abs(out[0][m][k]-out[1][m][k]) for m in out[1])
    if err>6e-13:
        raise ArithmeticError('The lower-cost comparison did not resolve between panel meshes')
    tag={} if observable=='forward_number' else dict(observable=observable)
    return dict(grid='fine',x=.05,charge_div=2,scalar_levels=levels,**tag,
                left_right_panel_edges=[e.tolist(),f.tolist()],**out[1],independent_panel_mesh_optimum_difference_q=err,
                scope='Conservative interval-support Wasserstein cost from a finite positive panel coupling. No intrapanel location, continuum optimum or physical momentum transfer inferred.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('left','right','states'):
        p.add_argument('--'+k,nargs=2,type=Path,required=True)
    p.add_argument('--levels',nargs=2,type=int,required=True)
    p.add_argument('--archive',type=Path,required=True)
    p.add_argument('--observable',choices=('forward_number','inverse_number','forward_excess'),default='forward_number')
    p.add_argument('--independent-bracket',action='store_true')
    a=p.parse_args()
    print(json.dumps(run((a.left,a.right),a.states,a.archive,a.levels,a.observable,independent_bracket=a.independent_bracket),indent=2,allow_nan=False))
