"""Construct small Pauli-safe sign-reversal witnesses in the fixed-moment aperture."""
import argparse
import json
from math import fsum

import numpy as np
from scipy.optimize import linprog
from scipy.special import expit
from scipy import sparse

from early_pole_daughter_frozen_moment_envelope import cell_objects
from early_pole_daughter_local_fraction_curvature_attribution import initial


def prepare(archive,sub):
    C,n,(y0,_,y1),_,depth=initial(archive)
    M=C.M
    e,f=M.weights(y0,True)[:2]
    a=np.asarray(M.edge_sets[2])
    active=np.flatnonzero((a[:-1]>=1.4138888888888888-2e-14)&(a[1:]<=1.4666666666666666+2e-14))
    sections=[]
    anchor_terms=[]
    k=C.scalar(n)
    for j in M.r[2].start+active:
        sample,A,W,total,_,_=cell_objects(C,int(j),e,f,sub,k,n,True)
        w,features,v,p=sample(24)
        lo,hi=sample.radial_bounds
        B=M.cells[j]
        fits=[B.fit(yy[j])[0] for yy in (y0,y1)]
        for sign,logits,yy in zip((1.,-1.),fits,(y0,y1),strict=True):
            calc=lambda q:expit(logits[0]+logits[1]*B.basis(q))
            F=calc(p)
            m=np.r_[yy[j],np.sum(w*features[2]*F)]
            ends=np.asarray([calc(lo),calc(hi)])
            fmin,fmax=ends.min(axis=0),ends.max(axis=0)
            sections.append(dict(sign=sign,A=A/total,m=m/total,signed_W=sign*W,fmin=fmin,fmax=fmax))
            anchor_terms.append(sign*float(np.sum(w*v*F)))
    return depth,len(n),sections,anchor_terms


def verify_repaired(z,sections):
    ends=np.r_[0,np.cumsum([len(v['signed_W']) for v in sections])]
    max_amp=0.
    defect=0.
    rates=[]
    for v,a,b in zip(sections,ends[:-1],ends[1:],strict=True):
        x=z[a:b].copy()
        A,m=v['A'],v['m']
        L=(A/A[0]).T
        lifted=L@np.linalg.inv(A@L)
        anchor=lifted@m
        if not np.all((anchor>0.)&(anchor<1.)):
            raise ValueError('Feasibility anchor is outside Pauli box.')
        for slack in (0.,1e-8,1e-7,1e-6,1e-5,1e-4,.001):
            y=(1.-slack)*x+slack*anchor
            y+=lifted@(m-A@y)
            if np.all((y>=0.)&(y<=1.)):
                break
        else:
            raise ValueError('Could not restore exact bounded moment witness.')
        defect=max(defect,float(max(abs(A@y-m))))
        rates.append(float(v['signed_W']@y))
        max_amp=max(max_amp,float(max(np.maximum(y-v['fmin'],v['fmax']-y))))
    if defect>5e-13:
        raise ValueError('Feasible witness moment residual too large.')
    return max_amp,defect,fsum(rates)


def optimize(sections,constant,target=1e-5):
    N=sum(len(v['signed_W']) for v in sections)
    init=np.concatenate([v['signed_W'] for v in sections])
    cap=max(abs(init).sum(),1e-8)
    high=np.concatenate([v['fmin'] for v in sections])
    low=np.concatenate([v['fmax'] for v in sections])
    one=sparse.csc_matrix(-np.ones((N,1)))
    mat=sparse.vstack((sparse.hstack((sparse.eye(N),one)),
                       sparse.hstack((-sparse.eye(N),one)),
                       sparse.csr_matrix(np.r_[-init/cap,0][None,:])),format='csc')
    b=np.r_[high,-low,(constant-target)/cap]
    eq=sparse.block_diag([v['A'] for v in sections],format='csc')
    E=sparse.hstack((eq,sparse.csc_matrix((eq.shape[0],1))),format='csc')
    r=linprog(np.r_[np.zeros(N),1.],A_ub=mat,b_ub=b,A_eq=E,
              b_eq=np.concatenate([v['m'] for v in sections]),bounds=(0.,1.),method='highs',
              options=dict(primal_feasibility_tolerance=1e-9,dual_feasibility_tolerance=1e-9))
    if r.status==2:
        return dict(status='infeasible_in_this_finite_approximation')
    if not r.success:
        raise ValueError(r.message)
    amp,defect,raw=verify_repaired(r.x[:-1],sections)
    p=constant+raw
    if p<=0.:
        raise ValueError('Repaired witness did not reverse the paired rate.')
    return dict(status='verified_positive_witness',solver_amplitude=float(r.x[-1]),
                verified_pointwise_amplitude=amp,positive_paired_rate=p,max_normalized_moment_residual=defect,
                selected_step_subcells=N,target_positive_rate=target)


def run(archive,sub):
    depth,n,sections,terms=prepare(archive,sub)
    native=fsum(x for x,v in zip(terms,sections,strict=True) if v['sign']==1.)
    local=fsum(x for x,v in zip(terms,sections,strict=True) if v['sign']==-1.)
    return dict(scope='Fixed frozen packet geometry and fixed number, energy, reconstructed reciprocal moment per daughter cell.',
                hard_depth=depth,cold_modes=n,subdivision=sub,fitted_pair=fsum(terms),
                native_fixed_influence=native,local_fixed_influence=-local,
                both_arms=optimize(sections,0.),
                native_only=optimize([x for x in sections if x['sign']==1.],local),
                local_only=optimize([x for x in sections if x['sign']==-1.],native))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('archives',nargs='+')
    p.add_argument('--subdivision',type=int,default=12)
    a=p.parse_args()
    if a.subdivision<1:
        p.error('Positive subdivision required.')
    print(json.dumps([run(f,a.subdivision) for f in a.archives],indent=2,allow_nan=False))
