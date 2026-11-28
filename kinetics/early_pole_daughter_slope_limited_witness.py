"""Continuous occupation witnesses on a uniform node grid, with explicit slope cap."""
import argparse
import json

import numpy as np
from scipy import sparse
from scipy.optimize import linprog
from scipy.special import expit

from early_pole_daughter_frozen_moment_envelope import cell_objects
from early_pole_daughter_local_fraction_curvature_attribution import initial


def fixture(archive,sub,fix_gauge=True):
    C,n,(y0,_,y1),_,depth=initial(archive)
    M=C.M
    e,f=M.weights(y0,True)[:2]
    pole=C.scalar(n)
    active=[j for j in range(M.r[2].start,M.r[2].stop)
            if M.cells[j].lo>=1.4138888888888888-2e-14 and M.cells[j].hi<=1.4666666666666666+2e-14]
    nvar=len(active)*sub+1
    K=np.zeros(nvar)
    dims=3 if fix_gauge else 2
    A=np.zeros((dims*len(active),nvar))
    m=np.zeros((2,len(A)))
    moment_interp=np.zeros(len(A))
    influence_interp=0.
    nodes=[]
    max_fitted_derivative=0.
    for k,j in enumerate(active):
        B=M.cells[j]
        h=np.linspace(B.lo,B.hi,sub+1)
        nodes.extend(h if k==0 else h[1:])
        sample,*_=cell_objects(C,j,e,f,1,pole,n,fix_gauge,h)
        w,features,influence,p=sample(24)
        ix=np.searchsorted(h,p)-1
        if np.any((ix<0)|(ix>=sub)):
            raise ValueError('Uncovered uniform shape node support.')
        frac=(p-h[ix])/(h[ix+1]-h[ix])
        def at_nodes(x):
            u=(w*x).ravel()
            return (np.bincount(ix.ravel(),weights=u*(1.-frac.ravel()),minlength=sub+1)
                    +np.bincount((ix+1).ravel(),weights=u*frac.ravel(),minlength=sub+1))
        rows=np.asarray([at_nodes(v) for v in features])
        size=rows[0].sum()
        ab=slice(k*sub,k*sub+sub+1)
        A[dims*k:dims*(k+1),ab]=rows/size
        K[ab]+=at_nodes(influence)
        # For any Lipschitz L occupation the chord interpolation
        # error is at most 2*L*h*t*(1-t).
        error1=2.*(h[ix+1]-h[ix])*frac*(1.-frac)
        moment_interp[dims*k:dims*(k+1)]=np.sum(
            abs(features)*(w*error1)[None,:,:],axis=(1,2))/size
        influence_interp+=float(np.sum(w*abs(influence)*error1))
        for i,y in enumerate((y0,y1)):
            z=B.fit(y[j])[0]
            F=expit(z[0]+z[1]*B.basis(p))
            derivative=np.abs(z[1]*p/np.hypot(p,B.m)/B.de*F*(1.-F))
            max_fitted_derivative=max(max_fitted_derivative,float(derivative.max()))
            target=np.r_[y[j],np.sum(w*features[2]*F)] if fix_gauge else y[j]
            m[i,dims*k:dims*(k+1)]=target/size
    exterior=[float(expit(f[active[0]-1,0]+f[active[0]-1,1])),
              float(expit(f[active[-1]+1,0]-f[active[-1]+1,1]))]
    return depth,len(n),A,m,K,np.asarray(nodes),max_fitted_derivative,moment_interp,influence_interp,exterior


def optimize(A,m,K,nodes,L,exterior=None):
    n=len(nodes)
    diff=np.diff(nodes)
    D=sparse.diags([-1./diff,1./diff],[0,1],shape=(n-1,n),format='csc')
    W=sparse.vstack((D,-D),format='csc')
    eye=sparse.eye(n)
    bd=sparse.vstack((eye,-eye,W),format='csc')
    free=np.arange(n) if exterior is None else np.arange(1,n-1)
    slack=np.ones(n)
    if exterior is not None:
        slack[[0,-1]]=0.
    bounds=[(0.,1.)]*n
    if exterior is not None:
        bounds[0]=(exterior[0],exterior[0])
        bounds[-1]=(exterior[1],exterior[1])
    slackscale=np.r_[slack,slack,np.full(2*(n-1),L)]
    tgt=np.r_[np.ones(n),np.zeros(n),np.full(2*(n-1),L)]
    X=sparse.hstack((bd,sparse.csc_matrix(slackscale[:,None])),format='csc')
    E=np.c_[A,np.zeros(len(A))]
    opt=dict(primal_feasibility_tolerance=1e-9,dual_feasibility_tolerance=1e-9)
    y=linprog(np.r_[np.zeros(n),-1.],A_ub=X,b_ub=tgt,A_eq=E,b_eq=m,
              bounds=bounds+[(0.,1.)],method='highs',options=opt)
    if y.status==2:
        return dict(status='no_feasible_profile_at_this_cap_on_this_grid')
    if not y.success:
        raise ValueError(y.message)
    A_free=A[:,free]
    G=np.zeros_like(A.T)
    G[free]=A_free.T@np.linalg.inv(A_free@A_free.T)
    anchor=y.x[:-1]
    anchor+=G@(m-A@anchor)
    scale=max(K.sum(),1e-10)
    ans=[]
    for s in (1.,-1.):
        r=linprog(s*K/scale,A_ub=W,b_ub=np.full(2*(n-1),L),A_eq=A,b_eq=m,
                  bounds=bounds,method='highs',options=opt)
        if not r.success:
            raise ValueError(r.message)
        for eps in (1e-7,1e-6,1e-5,1e-4,1e-3,.01):
            z=(1.-eps)*r.x+eps*anchor
            z+=G@(m-A@z)
            slope=float(max(abs(np.diff(z)/diff)))
            if np.all((z>=0.)&(z<=1.)) and slope<=L:
                break
        else:
            return dict(status='numerical_witness_repair_failed')
        defect=float(max(abs(A@z-m)))
        if defect>5e-13:
            raise ValueError('Numerical moment repair failed.')
        ans.append(dict(influence=float(K@z),numeric_finite_lp_optimum=float(s*r.fun*scale),verified_max_slope=slope,
                        normalized_moment_defect=defect))
    return dict(status='verified',minimum=ans[0],maximum=ans[1])


def relaxed_outer(A,m,K,nodes,L,moment_error,influence_error):
    step=np.diff(nodes)
    D=sparse.diags([-1./step,1./step],[0,1],shape=(len(step),len(nodes)),format='csc')
    M=sparse.vstack((A,-A,D,-D),format='csc')
    b=np.r_[m+L*moment_error,-m+L*moment_error,np.full(2*len(step),L)]
    scale=max(K.sum(),1e-9)
    v=[]
    for s in (1.,-1.):
        r=linprog(s*K/scale,A_ub=M,b_ub=b,bounds=(0.,1.),method='highs',
                  options=dict(primal_feasibility_tolerance=1e-9,dual_feasibility_tolerance=1e-9))
        if r.status==2:
            return None
        if not r.success:
            raise ValueError(r.message)
        multiplier=np.minimum(r.ineqlin.marginals,0.)
        leftover=s*K/scale-M.T@multiplier
        lower=float(b@multiplier+np.minimum(leftover,0.).sum())
        if lower>r.fun+1e-10:
            raise ValueError('Weak-dual check exceeds the relaxed finite optimum.')
        v.append(s*(lower*scale-L*influence_error))
    return v


def run(archive,sub,caps,fix_gauge=True,outer=False,match_exterior=False):
    d,q,A,m,K,nodes,fit,em,ew,edge=fixture(archive,sub,fix_gauge)
    if outer and match_exterior:
        raise ValueError('Use separate tests for the numerical outer relaxation and exterior matching.')
    out=[]
    for L in caps:
        x,y=(optimize(A,b,K,nodes,L,edge if match_exterior else None) for b in m)
        pair=None
        if x['status']==y['status']=='verified':
            pair=[x['minimum']['influence']-y['maximum']['influence'],
                  x['maximum']['influence']-y['minimum']['influence']]
        row=dict(slope_cap=L,native=x,local=y,paired_witness_range=pair)
        if outer:
            a,b=(relaxed_outer(A,z,K,nodes,L,em,ew) for z in m)
            row['numerical_interpolation_outer_range']=None if a is None or b is None else [a[0]-b[1],a[1]-b[0]]
            row['assumed_chord_influence_error_per_arm']=ew*L
            if pair is not None and row['numerical_interpolation_outer_range'] is not None:
                if any(a>b+1e-8 for a,b in zip((row['numerical_interpolation_outer_range'][0],pair[1]),
                                              (pair[0],row['numerical_interpolation_outer_range'][1]))):
                    raise ValueError('A verified piecewise witness fell outside the relaxed outer interval.')
        out.append(row)
    return dict(scope='Continuous Lipschitz-limited shapes on finite positive nodes in the original interval; not a proof of infeasibility outside this discretization.',
                hard_depth=d,cold_modes=q,subdivision=sub,fixed_gauge_moment=fix_gauge,
                sampled_max_within_cell_fitted_slope=fit,cases=out,
                **(dict(fixed_exterior_occupations=edge) if match_exterior else {}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('archives',nargs='+')
    p.add_argument('--subdivision',type=int,default=32)
    p.add_argument('--caps',type=float,nargs='+',default=[.5,1,2,4,8,16,32,64,128,256,512,1024])
    p.add_argument('--two-moments-only',action='store_true')
    p.add_argument('--with-outer',action='store_true')
    p.add_argument('--match-exterior',action='store_true')
    a=p.parse_args()
    if a.subdivision<4 or min(a.caps)<=0.:
        p.error('Require a positive slope and valid node intervals.')
    print(json.dumps([run(f,a.subdivision,a.caps,not a.two_moments_only,a.with_outer,a.match_exterior) for f in a.archives],indent=2,allow_nan=False))
