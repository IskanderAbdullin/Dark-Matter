"""Joint finite profile error when the absolute forward band totals are equal."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

from early_pole_endpoint_cold_joint_absolute_profile import make,gradients,solve as unconstrained
from early_pole_endpoint_cold_normalized_profile_envelope import fit_state,geometry


class AnchorUnavailable(ArithmeticError):
    def __init__(self, interval):
        self.interval=interval
        super().__init__('No strict numerical support for matching the forward totals')


def solve(a,b,field,member):
    """Match a common forward total without freely renormalizing other fields."""
    if field[0]!=0:
        raise ValueError('The anchored forward channel must lead the objective')
    n=a[0].shape[1]
    totals=a[1].sum(axis=1)
    ix=np.concatenate([np.arange(i*n,(i+1)*n) for i in field])
    if min(totals)<=0:
        raise ArithmeticError('Need a positive reference')
    scale=np.repeat(totals[list(field)],n)
    U,F,P=gradients(a,member,ix,scale)
    V,G,Q=gradients(b,member,ix,scale)
    D=U-V
    k,l=F.shape[1],G.shape[1]
    C=np.block([[P,np.zeros((len(P),l))],[np.zeros((len(Q),k)),Q]])
    J=np.c_[F,-G]
    w=J[:n].sum(axis=0)
    d=fsum(D[:n])
    g=np.array([d,d])
    corners=[]
    for i,side in enumerate((np.argmin,np.argmax)):
        r=np.zeros(k+l)
        for row in C:
            z=np.flatnonzero(row)
            r[z[side(w[z])]]=1.
        g[i]+=w@r
        corners.append(r)
    if not g[0]<-1e-10 or not g[1]>1e-10:
        raise AnchorUnavailable(g.tolist())
    m=len(D)
    A=np.block([[J,-np.eye(m)],[-J,-np.eye(m)]])
    bnd=np.r_[-D,D]
    E=np.block([[C,np.zeros((len(C),m))],[w[None,:],np.zeros((1,m))]])
    rhs=np.r_[np.ones(len(C)),-d]
    c=np.r_[np.zeros(k+l),np.ones(m)]
    res=linprog(c,A_ub=A,b_ub=bnd,A_eq=E,b_eq=rhs,bounds=(0.,None),method='highs',
        options={'primal_feasibility_tolerance':1e-10,'dual_feasibility_tolerance':1e-10})
    if not res.success:
        raise ArithmeticError(res.message)
    raw=res.x
    v=np.maximum(raw[:k+l],0.)
    t=C@v
    if min(t)<=0:
        raise ArithmeticError('Empty approximate feasible parent')
    v/=t@C
    r=d+w@v
    corner=int(r<0.)
    alpha=0. if r==0. else r/(r-g[corner])
    if not 0.<=alpha<1.:
        raise ArithmeticError('Invalid correction to the forward equality')
    z=(1-alpha)*v+alpha*corners[corner]
    L=U+F@z[:k]
    R=V+G@z[k:k+l]
    nd=L-R
    grouped=nd.reshape(len(field),n)
    loss=fsum(abs(nd))
    y0,y1=res.eqlin.marginals,res.ineqlin.marginals
    s=np.clip(y1[m:]-y1[:m],-1.,1.)
    s[:n]-=y0[-1]
    p=s@J
    support=fsum(s*D)+fsum(min(p[row==1.]) for row in C)
    residual=max(abs(fsum(nd[:n])),float(np.max(abs(C@z-1))),abs(loss-res.fun),
                 abs(loss-support),float(np.max(abs(E@raw-rhs))),
                 float(max(0.,np.max(A@raw-bnd),-np.min(raw))))
    if residual>2e-9 or support>loss+2e-12:
        raise ArithmeticError(('Cannot validate anchored optimum',residual))
    return dict(minimum_finite_loss=float(res.fun),projected_upper_witness=loss,
                support_lower_bound=support,unconstrained_forward_total_range=g.tolist(),
                absolute_witness_losses_by_field=np.sum(abs(grouped),axis=1).tolist(),
                signed_witness_differences_by_field=np.sum(grouped,axis=1).tolist(),
                matched_forward_totals_in_S2_scale=[fsum(L[:n]),fsum(R[:n])],
                maximum_defect=residual)


def run(left,right,states,archive):
    charts=[fit_state(p,archive,n,raw) for p,n,raw in zip(states,(2,4),(left,right),strict=True)]
    d,f,m,_=charts[0]
    lo,hi=1e-6,1.5e-5
    ce=np.r_[lo,d.ds.ed[(d.ds.ed>lo)&(d.ds.ed<hi)],hi]
    share=ce,f[np.searchsorted(d.ds.ed,(ce[:-1]+ce[1:])/2)-1],m
    out=[]
    for i in (0,1):
        Z=list(map(geometry,(left[i],right[i])))
        bins=np.intersect1d(Z[0][0],Z[1][0])
        if len(bins)<5 or [bins[0],bins[-1]]!=[lo,hi]:
            raise ArithmeticError('Need the exact held common partition')
        a,b=[make(*x,s,bins,share) for x,s in zip(Z,charts,strict=True)]
        rows={}
        for j,key in enumerate(('count_only','count_and_mean')):
            card={}
            for field,name in [((0,),'forward_only'),((0,1,2),'joint_three_fields')]:
                r=solve(a,b,field,j)
                u=unconstrained(a,b,field,j)
                if r['minimum_finite_loss']+2e-9<u['minimum_l1']:
                    raise ArithmeticError('Anchored optimum cannot improve relaxation')
                card[name]=dict(**r,without_total_anchor=u['minimum_l1'])
            rows[key]=card
        out.append((bins.tolist(),rows))
    (l,a),(r,b)=out
    err=max(abs(a[k][j]['minimum_finite_loss']-b[k][j]['minimum_finite_loss']) for k in a for j in a[k])
    if l!=r or err>2e-9:
        raise ArithmeticError('Anchor does not resolve across accepted meshes')
    return dict(x=.05,common_q_edges=r,results=b,max_paired_mesh_defect=err,
                scope='Finite common S2 cell moments, with matched forward band total. Chart optical fits and charge states retained. No continuum or physical inference.')


def refined_blends(left,right,baseline,states,archive,splits):
    from early_pole_endpoint_cold_shared_reference_scan import fit_refined_left,project,vol
    L=fit_refined_left(left,baseline,states[0],archive,splits)
    R=fit_state(states[1],archive,4,right)
    d,_,m,_=L
    lo,hi=1e-6,1.5e-5
    ce=np.r_[lo,d.ds.ed[(d.ds.ed>lo)&(d.ds.ed<hi)],hi]
    v=vol(ce[:-1],ce[1:])
    (c0,e0),(c1,e1)=map(lambda t:project(t,ce,m),(L,R))
    out=[]
    for i in (0,1):
        books=list(map(geometry,(left[i],right[i])))
        bins=np.intersect1d(books[0][0],books[1][0])
        if not np.array_equal(bins,np.intersect1d(geometry(baseline[i])[0],books[1][0])):
            raise ArithmeticError('Comparison bins changed under targeted reintegration')
        vals=[]
        for t in np.linspace(0.,1.,5):
            n=(1-t)*c0+t*c1
            e=(1-t)*e0+t*e1
            mean=np.divide(e,n,out=np.zeros_like(n),where=n>0)
            a,b=[make(*p,s,bins,(ce,n/v,m,mean)) for p,s in zip(books,(L,R),strict=True)]
            row=dict(S4_reference_weight=float(t))
            for field,name in [((0,),'forward_only'),((0,1,2),'joint_three_fields')]:
                base=unconstrained(a,b,field,1)['minimum_l1']
                try:
                    item=solve(a,b,field,1)
                    if item['minimum_finite_loss']<base-2e-9:
                        raise ArithmeticError('Exact match beat the relaxation')
                    row[name]=dict(without_total_anchor=base,**item)
                except AnchorUnavailable as err:
                    row[name]=dict(without_total_anchor=base,unconstrained_forward_total_range=err.interval,
                                   unresolved_or_unattainable=True)
            vals.append(row)
        out.append(vals)
    err=0.
    for a,b in zip(*out,strict=True):
        for k in ('forward_only','joint_three_fields'):
            x,y=a[k],b[k]
            if x.keys()!=y.keys():
                raise ArithmeticError('Refined meshes dispute attainability')
            err=max(err,np.max(abs(np.array(x['unconstrained_forward_total_range'])-y['unconstrained_forward_total_range'])))
            if 'minimum_finite_loss' in x:
                err=max(err,abs(x['minimum_finite_loss']-y['minimum_finite_loss']))
    if err>2e-9:
        raise ArithmeticError('Need pair-resolved forward overlap and joint optimum')
    return dict(x=.05,validation='Five predeclared mixtures on the full 17 parents after actual five-panel S2 reintegration',
                results=out[1],max_paired_mesh_defect=float(err),
                scope='Same held-response condition as native panel model; not a physical rate correction.')


if __name__=='__main__':
    q=argparse.ArgumentParser(description=__doc__)
    for name in ('left','right','states'):
        q.add_argument('--'+name,nargs=2,type=Path,required=True)
    q.add_argument('--archive',type=Path,required=True)
    q.add_argument('--left-baseline',nargs=2,type=Path)
    q.add_argument('--left-split-plan',type=Path)
    a=q.parse_args()
    if (a.left_baseline is None)!=(a.left_split_plan is None):
        q.error('For targeted reintegration include both baseline meshes and the split plan')
    if a.left_baseline is None:
        r=run(a.left,a.right,a.states,a.archive)
    else:
        splits=json.loads(a.left_split_plan.read_text())['targeted_parent_panel_splits']
        r=refined_blends(a.left,a.right,a.left_baseline,a.states,a.archive,splits)
    print(json.dumps(r,indent=2,allow_nan=False))
