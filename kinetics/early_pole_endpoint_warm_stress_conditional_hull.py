"""Count/energy-constrained pressure and pseudo-pressure envelopes at formal masses."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from early_pole_endpoint_gauss_warm_comparison import chart_for,gauss
from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack


def values(s,m):
    q=np.sqrt(s)
    E=np.hypot(q,m)
    v=q/E
    p=q*v/3
    return np.stack((p,p*v*v),axis=-1)


def bound(a,b,m,m0,n):
    """Global bound on each interpolation's absolute second-order remainder."""
    ea,eb=np.hypot([a,b],m0)
    h=(b-a)*(b+a)/(eb+ea)
    if m==0:
        return np.full(2,h*h*m0*m0/(24*a**3*n*n))
    if m==m0:
        c=m0*m0
        return h*h/(24*n*n)*np.asarray([2*c/ea**3,4*c/ea**3+12*c*c/ea**5])
    lo,hi=a*a+m*m,b*b
    out=[]
    for k in (1,2):
        p=k-.5
        x=k*hi**(k-1)/lo**p+p*hi**k/lo**(p+1)
        y=(0. if k==1 else k*(k-1)*hi**(k-2)/lo**p)
        y+=2*k*p*hi**(k-1)/lo**(p+1)+p*(p+1)*hi**k/lo**(p+2)
        # One third is part of both pressure kernels.
        out.append(h*h*(4*eb*eb*y+2*x)/(24*n*n))
    return np.asarray(out)


def hull(a,b,m,m0,mean,n):
    if not (0<a<b and 0.<=mean<=1. and n>=2):
        raise ArithmeticError('Invalid positive retained conditional energy')
    u=np.linspace(0.,1.,n+1)
    ea,eb=np.hypot([a,b],m0)
    h=(b-a)*(b+a)/(eb+ea)
    s=a*a+(u*h)*(2*ea+u*h)
    M=values(s,m)
    i=np.flatnonzero(u<mean); j=np.flatnonzero(u>mean)
    choices=[]
    if len(i) and len(j):
        i=i[:,None]; j=j[None,:]
        t=((mean-u[i])/(u[j]-u[i]))[...,None]
        X=(1-t)*M[i]+t*M[j]
        choices.append(X.reshape(-1,2))
    if np.any(u==mean):
        choices.append(M[u==mean])
    Z=np.concatenate(choices)
    l=np.min(Z,axis=0)
    r=np.max(Z,axis=0)
    e=bound(a,b,m,m0,n)
    return np.maximum(0,l-e),r+e,e,np.stack((l,r))


def conditioned(base,card,M0):
    gauss(card,base,M0,.05) # strict state and observation checks
    tag,d,repeats=chart_for(card)
    e=d.ds.ed
    y=np.asarray(card['nodal_occupation'])
    w=d.ds.w*y[:d.ns]
    k=np.r_[0,np.cumsum(repeats)]
    output=[]
    for i,(a,b) in enumerate(zip(e[:-1],e[1:],strict=True)):
        j=slice(k[i],k[i+1]); q=d.ds.p[j]
        N=fsum(w[j])
        ha=np.hypot(a,M0)
        E=(b-a)*(b+a)/(np.hypot(b,M0)+ha)
        u=(q-a)*(q+a)/(np.hypot(q,M0)+ha)/E
        mu=0. if N==0 else fsum(w[j]*u)/N
        if not 0.<=mu<=1.:
            raise ArithmeticError('Nodal conditional energy is outside its own cell')
        output.append((a,b,N,mu,q,w[j],repeats[i]>1))
    return tag,output


def measure(c,m,m0,n):
    G,W=roots_legendre(96)
    bounds=[]; witness=[]; errors=[]; feasible=[]
    for a,b,N,mean,q,w,warm in c:
        if warm:
            lo,hi,e,inner=hull(a,b,m,m0,mean,n)
            p=w@values(q*q,m)
        else:
            v=a+(b-a)*(G+1)/2
            V=W*v*v
            p=N*(V@values(v*v,m)/fsum(V))
            lo=hi=np.zeros(2) if N==0 else p/N
            e=np.zeros(2)
            inner=np.stack((lo,hi))
        if np.any(p<N*lo-5e-14*np.maximum(p,1e-99)) or np.any(p>N*hi+5e-14*np.maximum(p,1e-99)):
            raise ArithmeticError('Direct nodal positive witness escaped the conditional envelope')
        bounds.append(np.stack((N*lo,N*hi)))
        witness.append(p)
        errors.append(N*e)
        feasible.append(N*inner)
    return tuple(np.asarray(x).sum(axis=0) for x in (bounds,witness,errors,feasible))


def run(base,states,receipts,levels=(3,4)):
    levels=tuple(levels)
    if levels not in ((3,4),(4,5),(5,6)):
        raise ArithmeticError('Need two independent adjacent accepted warm Gauss histories')
    obs=f'recorded_G{levels[1]}_G{levels[0]}_ratio_minus_one'
    d,y=unpack(base,False,.05)
    m0=d.state(.05,y)[0]
    C=[conditioned(d,g,m0) for g in states]
    if tuple(a for a,_ in C)!=levels:
        raise ArithmeticError('Wrong ordered positive warm archives')
    check=[[s['rows'][-1]['totals'][2:] for s in r['measurements']] for r in receipts]
    mm=[[s['test_mass'] for s in r['measurements']] for r in receipts]
    if mm[0]!=mm[1] or any(r['reference_fitted_mass']!=m0 for r in receipts):
        raise ArithmeticError('Frozen comparison receipt mismatch')
    out=[]
    for j,m in enumerate(mm[0]):
        u,v=[[measure(c,m,m0,n) for _,c in C] for n in (64,256)]
        delta=[]
        for k,((lo,hi),w,e,inner) in enumerate(v):
            x=check[k][j]
            if np.any(np.asarray(x)<lo-2e-14*lo) or np.any(np.asarray(x)>hi+2e-14*hi):
                raise ArithmeticError('The independently archived exponential witness escaped')
            a=u[k][0]
            if np.any(lo>a[1]) or np.any(hi<a[0]):
                raise ArithmeticError('Independent hull meshes do not overlap')
            delta.append(dict(chart=C[k][0],conditional_interval=[lo.tolist(),hi.tolist()],
                              finite_atom_extrema=inner.tolist(),interpolation_remainder=e.tolist(),direct_witness=w.tolist(),reconstructed_witness=x))
        a,b=[p[0] for p in v]
        ratio=np.array([b[0]/a[1]-1,b[1]/a[0]-1]).T
        low,high=[p[3] for p in v]
        feasible=np.array([high[0]/low[1]-1,high[1]/low[0]-1]).T
        if np.any(feasible[:,0]<ratio[:,0]-2e-15) or np.any(feasible[:,1]>ratio[:,1]+2e-15):
            raise ArithmeticError('Constructive atom witness escaped its outer envelope')
        out.append(dict(formal_mass=m,**{obs:(np.asarray(check[1][j])/check[0][j]-1).tolist()},
                        conditional_ratio_minus_one=ratio.tolist(),feasible_atom_ratio_minus_one=feasible.tolist(),charts=delta))
    return dict(x=.05,reference_mass=m0,field_order=['pressure','pseudo_pressure'],rows=out,
                scope='All nonnegative in-cell measures preserving warm count and mean excess energy at the reference mass, cold kept in archived constant convention. Analytic interpolation remainder, ordinary floating point; no physical clock or continuum dynamics claim.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--native',type=Path,required=True)
    p.add_argument('--gauss',nargs=2,type=Path,required=True)
    p.add_argument('--receipts',nargs=2,type=Path,required=True)
    p.add_argument('--levels',type=int,nargs=2,default=(3,4))
    a=p.parse_args()
    print(json.dumps(run(read(a.native),[read(p) for p in a.gauss],[read(p) for p in a.receipts],a.levels),indent=2,allow_nan=False))
