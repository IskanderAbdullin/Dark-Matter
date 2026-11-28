"""Quantify departures from decreasing occupation without prescribing a slope."""
import argparse
import json

import numpy as np
from scipy import sparse
from scipy.optimize import linprog

from early_pole_daughter_slope_limited_witness import fixture


def solve(A,m,K,x,edge,target):
    n=len(x)
    step=sparse.diags([-1.,1.],[0,1],shape=(n-1,n),format='csc')
    D=sparse.block_diag((step,step),format='csc')
    d=D.shape[0]
    W=sparse.hstack((D,-sparse.eye(d)),format='csc')
    bounds=np.zeros(d)
    E0=sparse.block_diag((A,A),format='csc')
    E=sparse.hstack((E0,sparse.csc_matrix((E0.shape[0],d))),format='csc')
    if target is not None:
        z=max(K.sum(),1e-9)
        W=sparse.vstack((W,sparse.csr_matrix(np.r_[-K/z,K/z,np.zeros(d)][None,:])),format='csc')
        bounds=np.r_[bounds,-target/z]
    box=[(0.,1.)]*n
    if edge:
        box[0]=(edge[0],edge[0])
        box[-1]=(edge[1],edge[1])
    obj=np.r_[np.zeros(2*n),np.ones(d)]
    r=linprog(obj,A_ub=W,b_ub=bounds,A_eq=E,b_eq=m.ravel(),bounds=box*2+[(0.,1.)]*d,method='highs',
              options=dict(primal_feasibility_tolerance=1e-9,dual_feasibility_tolerance=1e-9))
    if not r.success:
        raise ValueError(r.message)
    lam=r.eqlin.marginals
    mu=np.minimum(r.ineqlin.marginals,0.)
    a=obj-E.T@lam-W.T@mu
    lower=float(lam@m.ravel()+mu@bounds+sum(z*(lo if z>=0 else hi)
                for z,(lo,hi) in zip(a,box*2+[(0.,1.)]*d,strict=True)))
    free=np.arange(n) if not edge else np.arange(1,n-1)
    Q=A[:,free]
    adj=np.zeros_like(A.T)
    adj[free]=Q.T@np.linalg.inv(Q@Q.T)
    repaired=[]
    for value,b in zip(r.x[:2*n].reshape(2,n),m,strict=True):
        v=value+adj@(b-A@value)
        if not np.all((v>=0.)&(v<=1.)):
            adj2=np.isin(np.arange(n),free)
            B=sparse.vstack((sparse.hstack((sparse.eye(n),sparse.csc_matrix(adj2[:,None]))),
                             sparse.hstack((-sparse.eye(n),sparse.csc_matrix(adj2[:,None])))),format='csc')
            w=linprog(np.r_[np.zeros(n),-1.],A_ub=B,b_ub=np.r_[np.ones(n),np.zeros(n)],
                      A_eq=np.c_[A,np.zeros(len(A))],b_eq=b,bounds=box+[(0.,1.)],method='highs')
            if not w.success:
                raise ValueError(w.message)
            anchor=w.x[:-1]
            for t in (1e-8,1e-7,1e-6,1e-5,1e-4):
                v=(1.-t)*value+t*anchor
                v+=adj@(b-A@v)
                if np.all((v>=0.)&(v<=1.)):
                    break
            else:
                raise ValueError('Moment repair left the Pauli box.')
        if max(abs(A@v-b))>5e-13 or (edge and max(abs(v[[0,-1]]-edge))>2e-15):
            raise ValueError('Bad repaired witness constraints.')
        repaired.append(v)
    pos=[float(np.maximum(np.diff(y),0.).sum()) for y in repaired]
    pair=float(K@(repaired[0]-repaired[1]))
    if sum(pos)<lower-1e-8 or (target is not None and pair<target/2):
        raise ValueError('Witness does not verify its reported variation or sign.')
    return dict(total_upward_variation=sum(pos),upward_variation_by_arm=pos,
                numerical_dual_lower=lower,verified_paired_rate=pair)


def run(archive,sub,gauge,ext):
    d,q,A,m,K,x,F,_,_,edge=fixture(archive,sub,gauge)
    edge=edge if ext else None
    return dict(scope='Only the fixed interval, exact modeled moments, held affine packet geometry.',
                hard_depth=d,cold_modes=q,subdivision=sub,gauge_moment=gauge,matched_exterior=ext,
                moment_feasibility=solve(A,m,K,x,edge,None),
                positive_witness=solve(A,m,K,x,edge,1e-5))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('archives',nargs='+')
    p.add_argument('--subdivision',type=int,default=64)
    p.add_argument('--match-exterior',action='store_true')
    p.add_argument('--two-moments-only',action='store_true')
    a=p.parse_args()
    print(json.dumps([run(f,a.subdivision,not a.two_moments_only,a.match_exterior) for f in a.archives],
                     indent=2,allow_nan=False))
