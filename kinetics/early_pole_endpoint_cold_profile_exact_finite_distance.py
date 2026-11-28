"""Minimize normalized profile separation with a jointly constrained finite linear program."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

from early_pole_endpoint_cold_normalized_profile_envelope import build,fit_state,geometry,compare


def as_system(card,mode):
    """Charnes-Cooper variables: common scale and shell-vertex weights."""
    V=card[mode]
    n=1+sum(len(z) for z in V)
    qs=np.zeros((len(card['base']),n))
    qs[:,0]=card['offset']/card['total']
    eq=np.zeros((len(V)+1,n))
    k=1
    for j,z in enumerate(V):
        l=k+len(z)
        eq[j,0]=-1.
        eq[j,k:l]=1.
        qs[:,k:l]=z.T/card['total']
        k=l
    eq[-1]=qs.sum(axis=0)
    target=np.r_[np.zeros(len(V)),1.]
    if np.min(qs)<0 or not qs[:,0].sum()+fsum(min(v.sum(axis=1))/card['total'] for v in V)>0:
        raise ArithmeticError('Need a strictly positive total forward coefficient')
    return qs,eq,target


def profile_from_projected_blocks(card,mode,proposal):
    """Normalize the response of independently projected positive simplex blocks."""
    v=np.asarray(card['offset'],dtype=float).copy()/card['total']
    k=1
    for z in card[mode]:
        a=np.maximum(0.,proposal[k:k+len(z)])
        N=fsum(a)
        if np.isfinite(N) and N>0.:
            a=a/N
        else:  # A valid witness even when the optimizer returns an empty block.
            a=np.full(len(z),1./len(z))
        v+=a@(z/card['total'])
        k+=len(z)
    if k!=len(proposal) or min(v)<0. or not np.all(np.isfinite(v)) or not fsum(v)>0.:
        raise ArithmeticError('Invalid projected normalized profile')
    return v/fsum(v)


def normalized_support_lower(card,mode,direction):
    """Independent lower bound on min(direction dot normalized response).

    With bin-response numerator U and total denominator T, the minimum
    of U - t*T separates over the positive simplex cells. For a trial t
    with a negative residual h, min(U/T) >= t + h/min(T).
    """
    V=[np.asarray(card['offset'],dtype=float)[None,:]/card['total']]
    V.extend(z/card['total'] for z in card[mode])
    blocks=[(v@direction,v.sum(axis=1)) for v in V]
    minimum=fsum(min(t) for _,t in blocks)
    maximum=fsum(max(t) for _,t in blocks)
    if not 0.<minimum<=maximum or not np.isfinite(maximum):
        raise ArithmeticError('The support probe needs a positive bounded response')
    lo,hi=float(min(direction)),float(max(direction))
    # Bisection only proposes a threshold. The separately evaluated
    # residual correction below determines the reported lower bound.
    for _ in range(50):
        pivot=(lo+hi)/2
        if fsum(min(u-pivot*t) for u,t in blocks)<0.:
            hi=pivot
        else:
            lo=pivot
    pivot=(lo+hi)/2
    h=fsum(min(u-pivot*t) for u,t in blocks)
    # A small scale-aware floating-point guard, not directed rounding.
    guard=64*np.finfo(float).eps*maximum
    return max(float(min(direction)),pivot+min(0.,h-guard)/minimum)


def verify_normalized_tv_bracket(a,b,mode,pa,pb,sign):
    sign=np.clip(sign,-.5,.5)
    lower=max(0.,normalized_support_lower(a,mode,sign)+normalized_support_lower(b,mode,-sign))
    a_profile=profile_from_projected_blocks(a,mode,pa)
    b_profile=profile_from_projected_blocks(b,mode,pb)
    upper=fsum(abs(a_profile-b_profile))/2+64*np.finfo(float).eps
    if not np.isfinite(upper) or upper<lower or upper-lower>2e-9:
        raise ArithmeticError(('Independent support / projected TV bracket did not close',lower,upper))
    return dict(support_lower_bound=lower,projected_profile_upper_bound=upper,
                verified_TV_bracket_width=upper-lower,signed_support=sign.tolist(),
                projected_profiles=[a_profile.tolist(),b_profile.tolist()])


def solve(a,b,mode,*,independent_bracket=False):
    F,P,x=as_system(a,mode)
    G,Q,y=as_system(b,mode)
    n,k=F.shape
    if not len(G)==n:
        raise ValueError('Need the exact same profile bins')
    j=G.shape[1]
    delta=np.c_[F,-G]
    A=np.block([[delta,-np.eye(n)],[-delta,-np.eye(n)]])
    E=np.block([[P,np.zeros((len(P),j+n))],
                [np.zeros((len(Q),k)),Q,np.zeros((len(Q),n))]])
    rhs=np.r_[x,y]
    c=np.r_[np.zeros(k+j),np.full(n,.5)]
    z=linprog(c,A_ub=A,b_ub=np.zeros(2*n),A_eq=E,b_eq=rhs,bounds=(0.,None),
              method='highs',options={'primal_feasibility_tolerance':1e-10,'dual_feasibility_tolerance':1e-10})
    if not z.success:
        raise ArithmeticError(('Finite normalized profile LP failed',z.message))
    w=z.x
    u,v=F@w[:k],G@w[k:k+j]
    direct=fsum(abs(u-v))/2
    eq=float(np.max(abs(E@w-rhs)))
    feas=max(float(np.max(A@w)),float(-min(w)),0.)
    yy=z.eqlin.marginals
    zz=z.ineqlin.marginals
    dual=max(float(np.max(E.T@yy+A.T@zz-c)),float(np.max(zz)),0.)
    lower=float(rhs@yy)
    gap=max(abs(float(z.fun)-lower),abs(direct-z.fun))
    if independent_bracket:
        check=verify_normalized_tv_bracket(a,b,mode,w[:k],w[k:k+j],zz[n:]-zz[:n])
        return dict(minimum_finite_common_bin_TV=(check['support_lower_bound']+check['projected_profile_upper_bound'])/2,
                    direct_witness_TV=check['projected_profile_upper_bound'],independent_profile_bracket=check,
                    optimizer_nominal_maximum_residual=max(eq,feas,dual,gap))
    if max(eq,feas,dual,gap)>2e-9 or min(u)<-2e-12 or min(v)<-2e-12:
        raise ArithmeticError('Primal/dual or independent profile witness failed')
    return dict(minimum_finite_common_bin_TV=float(z.fun),direct_witness_TV=direct,
                verified_dual_value=lower,primal_dual_and_witness_maximum_difference=gap,
                maximum_equality_defect=eq,maximum_inequality_or_positivity_defect=feas,
                maximum_dual_inequality_defect=dual,optimized_profiles=[u.tolist(),v.tolist()],
                chart_total_response_relative_to_separate_original=[1/w[0],1/w[k]])


def run(left,right,states,archive,levels,observable='forward_number',*,independent_bracket=False):
    if tuple(levels) not in ((2,4),(4,8),(8,16)):
        raise ArithmeticError('The accepted pair must be independent and adjacent')
    entries=[fit_state(p,archive,n,pair) for p,n,pair in zip(states,levels,(left,right),strict=True)]
    outputs=[]
    for m in (0,1):
        a,b=map(geometry,(left[m],right[m]))
        common=np.intersect1d(a[0],b[0])
        if len(common)<5 or [common[0],common[-1]]!=[1e-6,1.5e-5]:
            raise ArithmeticError('Need complete exactly shared cold band edges')
        p,q=[build(d,e,common,raw=True,observable=observable) for d,e in zip((a,b),entries,strict=True)]
        x,y=[build(d,e,common,observable=observable) for d,e in zip((a,b),entries,strict=True)]
        outer=compare(x,y,common)
        r={}
        for mode in ('count_only','count_and_mean_excess'):
            ans=solve(p,q,mode,independent_bracket=independent_bracket)
            lower=outer[mode]['robust_total_variation_lower_bound']
            base=fsum(abs(np.asarray(p['base'])-q['base']))/2
            if not lower-2e-9<=ans['minimum_finite_common_bin_TV']<=base+2e-9:
                raise ArithmeticError('Joint optimization not consistent with lower bound and base witness')
            r[mode]=dict(**ans,per_group_conservative_lower_bound=lower,original_step_TV=base)
        outputs.append((common.tolist(),r))
    (u,a),(v,b)=outputs
    err=max(abs(a[m]['minimum_finite_common_bin_TV']-b[m]['minimum_finite_common_bin_TV']) for m in b)
    if u!=v or err>2e-9:
        raise ArithmeticError('Finite optima disagree across the two accepted panel meshes')
    tag={} if observable=='forward_number' else dict(observable=observable)
    return dict(grid='fine',x=.05,charge_div=2,scalar_levels=levels,**tag,common_q_edges=v,
                **b,paired_mesh_absolute_optimum_defect=err,
                scope='Numerically resolved minimum on finite, conditional panel variables and exact common coarse q bins. Not the physical or continuum profile distance; no directed rounding.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('left','right','states'):
        p.add_argument('--'+k,nargs=2,type=Path,required=True)
    p.add_argument('--levels',nargs=2,type=int,required=True)
    p.add_argument('--archive',type=Path,required=True)
    p.add_argument('--observable',choices=('forward_number','inverse_number','forward_excess'),default='forward_number')
    p.add_argument('--independent-bracket',action='store_true')
    a=p.parse_args()
    print(json.dumps(run(a.left,a.right,a.states,a.archive,a.levels,a.observable,independent_bracket=a.independent_bracket),indent=2,allow_nan=False))
