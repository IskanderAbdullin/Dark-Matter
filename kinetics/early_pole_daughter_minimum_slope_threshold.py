"""Separate feasible moment interpolation from the slope needed for a positive witness."""
import argparse
import json

import numpy as np
from scipy import sparse
from scipy.optimize import linprog

from early_pole_daughter_slope_limited_witness import fixture
from early_pole_daughter_local_fraction_curvature_attribution import initial


def solve_joint(A,m,K,node,edge,target):
    N=len(node)
    dx=np.diff(node)
    D=sparse.diags([-1.,1.],[0,1],shape=(N-1,N),format='csc')
    laps=sparse.vstack((D,-D),format='csc')
    block=sparse.block_diag((laps,laps),format='csc')
    col=-512*np.tile(dx,4)
    limit=sparse.hstack((block,sparse.csc_matrix(col[:,None])),format='csc')
    b=np.zeros(block.shape[0])
    sc=max(K.sum(),1e-8)
    if target is not None:
        limit=sparse.vstack((limit,sparse.csr_matrix(np.r_[-K/sc,K/sc,0.][None,:])),format='csc')
        b=np.r_[b,-target/sc]
    abe=sparse.block_diag((A,A),format='csc')
    Ae=sparse.hstack((abe,sparse.csc_matrix((abe.shape[0],1))),format='csc')
    ib=[(0.,1.)]*N
    if edge:
        ib[0]=(edge[0],edge[0])
        ib[-1]=(edge[1],edge[1])
    bounds=ib+ib+[(0.,None)]
    opts=dict(primal_feasibility_tolerance=1e-9,dual_feasibility_tolerance=1e-9)
    obj=np.r_[np.zeros(2*N),512.]
    sol=linprog(obj,A_ub=limit,b_ub=b,A_eq=Ae,b_eq=m.ravel(),bounds=bounds,method='highs',options=opts)
    if not sol.success:
        raise ValueError(sol.message)
    # Evaluate a weak-dual lower bound independently of the return status.
    # Clamp inequality multipliers to the dual sign, then take the infimum
    # of the residual objective independently over each variable interval.
    lam=sol.eqlin.marginals
    mu=np.minimum(sol.ineqlin.marginals,0.)
    r=obj-Ae.T@lam-limit.T@mu
    dual=float(lam@m.ravel()+mu@b)
    # A finite cap upper bound can be set larger than any witness without
    # changing the minimum; evaluate its contribution with that bound.
    final_cap_upper=max(sol.x[-1]*2.,8.)
    for coef,(a,z) in zip(r,bounds[:-1]+[(0.,final_cap_upper)],strict=True):
        dual+=coef*(a if coef>=0. else z)
    f=np.arange(N) if not edge else np.arange(1,N-1)
    A0=A[:,f]
    lift=np.zeros_like(A.T)
    lift[f]=A0.T@np.linalg.inv(A0@A0.T)
    y=sol.x[:-1].reshape(2,N)
    repaired=[]
    defects=[]
    for x,goal in zip(y,m,strict=True):
        # Move to a small interior anchor only when necessary.
        yc=x+lift@(goal-A@x)
        if not np.all((yc>=0.)&(yc<=1.)):
            z=linprog(np.r_[np.zeros(N),-1.],A_eq=np.c_[A,np.zeros(len(A))],b_eq=goal,
                      A_ub=sparse.hstack((sparse.vstack((sparse.eye(N),-sparse.eye(N))),
                                         sparse.csc_matrix(np.tile(np.isin(np.arange(N),f),2)[:,None]))),
                      b_ub=np.r_[np.ones(N),np.zeros(N)],bounds=ib+[(0.,1.)],method='highs',options=opts)
            if not z.success:
                raise ValueError('Unable to find strict interior repair state.')
            anchor=z.x[:-1]
            for eps in (1e-8,1e-7,1e-6,1e-5,1e-4):
                yc=(1.-eps)*x+eps*anchor
                yc+=lift@(goal-A@yc)
                if np.all((yc>=0.)&(yc<=1.)):
                    break
            else:
                raise ValueError('No exact positive repair of minimum-slope nodal state.')
        if edge and np.max(abs(yc[[0,-1]]-edge))>2e-15:
            raise ValueError('Incorrect externally matched boundaries.')
        repaired.append(yc)
        defects.append(float(max(abs(A@yc-goal))))
    a,b=repaired
    paired=float(K@(a-b))
    vmax=float(max(max(abs(np.diff(v)/dx)) for v in (a,b)))
    if (max(defects)>5e-13 or (target is not None and paired<target/2)
            or dual>vmax+1e-7 or vmax>512.*final_cap_upper):
        raise ValueError('Minimum-slope candidate did not pass the independent verification.')
    return dict(numerical_finite_lower_slope=dual,verified_witness_slope=vmax,verified_paired_rate=paired,
                requested_positive_margin=target,max_normalized_moment_defect=max(defects))


def run(archive,sub,gauge,exterior):
    d,q,A,m,K,node,fit,em,ew,edges=fixture(archive,sub,gauge)
    ed=edges if exterior else None
    C,_,(native,_,local),_,_=initial(archive)
    M=C.M
    cells=[k for k in range(M.r[2].start,M.r[2].stop)
           if M.cells[k].lo>=1.4138888888888888-2e-14 and M.cells[k].hi<=1.4666666666666666+2e-14]
    low,high=(np.array([getattr(M.cells[k],v) for k in cells]) for v in ('lo','hi'))
    pbar=.75*(high**4-low**4)/(high**3-low**3)
    measure=2*(high**3-low**3)/(3*np.pi**2)
    floors=[]
    for y in (native,local):
        avg=y[cells,0]/measure
        F=list(abs(np.diff(avg))/np.diff(pbar))
        if ed is not None:
            F.extend(abs(avg-ed[0])/(pbar-low[0]))
            F.extend(abs(avg-ed[1])/(high[-1]-pbar))
        floors.append(float(max(F)))
    return dict(analytic_average_based_feasibility_lower_slope=floors,hard_depth=d,cold_modes=q,subdivision=sub,fixed_gauge_moment=gauge,match_exterior=exterior,
                minimum_feasibility=solve_joint(A,m,K,node,ed,None),
                minimum_positive=solve_joint(A,m,K,node,ed,1e-5))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('archives',nargs='+')
    p.add_argument('--subdivision',type=int,default=128)
    p.add_argument('--two-moments-only',action='store_true')
    p.add_argument('--match-exterior',action='store_true')
    a=p.parse_args()
    print(json.dumps([run(f,a.subdivision,not a.two_moments_only,a.match_exterior)
                      for f in a.archives],indent=2,allow_nan=False))
