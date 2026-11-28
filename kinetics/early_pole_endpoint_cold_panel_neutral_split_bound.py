"""Bound the correlated energy/number neutrality split on a finite cold paneled chart."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_cold_panel_reweight_envelope import mean_vertices,panel_excess
from early_pole_endpoint_cold_scalar_band_ladder import verify_pair
from early_pole_endpoint_continuum_cold_band_receipt import read as panels
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart
from early_pole_endpoint_stencil_remap_pilot import SPECS,load


def product_bounds(a,b):
    c=np.outer(a,b)
    return float(np.min(c)),float(np.max(c))


def bounds(offset,choices,alpha,beta,m):
    # Dimensionless vector b,d,x,y = B_N,D_N,B_X,D_X, each relative to its baseline.
    lo,hi=[offset+np.array([fsum(h(r[:,j]) for r in choices) for j in range(4)]) for h in (min,max)]
    if min(lo)<=0 or not np.all(np.isfinite([lo,hi])) or np.any(lo>1+4e-14) or np.any(hi<1-4e-14):
        raise ArithmeticError('Baseline left the feasible positive moments')
    den=(lo[0]*(m*lo[1]+beta*lo[3]),hi[0]*(m*hi[1]+beta*hi[3]))
    u=np.stack((lo,hi))-1
    p=product_bounds(u[:,2],u[:,1])
    q=product_bounds(u[:,0],u[:,3])
    # Expanding around the original four normalized moments preserves the
    # cancellation between the two bulk products up to the verified bounds.
    L=np.array([-beta,alpha,alpha,-beta])
    def linear(v):
        a=np.asarray(v)
        return [fsum([float((offset-1)@a),*(f(r@a) for r in choices)]) for f in (min,max)]
    lin=linear(L)
    uncoupled=(alpha*p[0]-beta*q[1],alpha*p[1]-beta*q[0])
    r0=alpha-beta
    # Rewrite the same quadratic exactly in terms of two correlated small
    # differences, rather than subtracting two independent bulk products.
    dx_y=linear([0,0,1,-1])
    dd_b=linear([-1,1,0,0])
    u1=product_bounds(u[:,1],dx_y)
    u2=product_bounds(u[:,3],dd_b)
    u3=[min(r0*q[0],r0*q[1]),max(r0*q[0],r0*q[1])]
    tight=(alpha*(u1[0]+u2[0])+u3[0],alpha*(u1[1]+u2[1])+u3[1])
    bil=(max(uncoupled[0],tight[0]),min(uncoupled[1],tight[1]))
    if bil[0]>bil[1]:
        raise ArithmeticError('Two independently conservative quadratic remainders did not intersect')
    src=[r0+u+v for u,v in zip(lin,bil,strict=True)]
    frac=np.asarray(src)[:,None]/np.asarray(den)[None,:]
    return dict(conditional_relative_split_interval=[float(np.min(frac)),float(np.max(frac))],
                linear_numerator_change_interval=lin,quadratic_numerator_change_interval=bil,
                independently_boxed_quadratic_interval=uncoupled,
                correlated_linear_deviations=[dx_y,dd_b],
                positive_denominator_interval=den,
                baseline_relative_split=r0/(m+beta))


def run(raw,state,archive,n):
    r=verify_pair(raw,state,archive,n)['diagnostic']
    entry,_=read(state)
    d=ScalarSubshellChart(load('fine')[0],SPECS['fine'],2,n,1e-5)
    m=r['full_pole_and_snapshot']['scalar_mass']
    def one(p):
        records=panels(p)
        e=np.array([records[0]['bounds'][0]]+[x['bounds'][1] for x in records])
        a,b=e[:-1],e[1:]
        g=np.searchsorted(d.ds.ed,(a+b)/2)-1
        if np.any((g<0)|(g>=d.ns)) or np.any((a<d.ds.ed[g])|(b>d.ds.ed[g+1])):
            raise ArithmeticError('An accepted q panel straddles a scalar shell')
        f=np.asarray(entry['endpoint_occupation'])[g]
        z=np.array([q['coeff'] for q in records])[:,[0,2],:]
        A=z[:,:,0]/(1+f[:,None])
        B=np.divide(z[:,:,1],f[:,None],out=np.zeros_like(A),where=f[:,None]>0)
        if min(A.ravel())<0 or np.any(B>A) or np.any((f==0)&np.any(z[:,:,1]!=0,axis=1)):
            raise ArithmeticError('The local affine panel components are not positive')
        D=A-B
        V=(b-a)*(b*b+b*a+a*a)/3
        X=panel_excess(a,b,m,48)
        col=np.stack((B[:,0],D[:,0],B[:,1],D[:,1]),axis=-1)/V[:,None]
        c=np.array([0,fsum(A[:,0]),0,fsum(A[:,1])])
        base=np.array([fsum(z[:,0,1]),fsum(z[:,0,2]),fsum(z[:,1,1]),fsum(z[:,1,2])])
        if min(base)<=0 or np.max(abs((c+(V*f)@col)/base-1))>5e-14:
            raise ArithmeticError('Finite-panel normalization does not close')
        U,W=[],[]
        for j in np.unique(g):
            i=g==j
            w=V[i]*f[i]
            N=fsum(w)
            if N==0:
                continue
            q=col[i]/base
            t=fsum(w*X[i])/N
            mat=mean_vertices(X[i],t) if len(w)>1 else np.ones((1,1))
            U.append(N*q)
            W.append(N*(mat@q))
        p,q=[bounds(c/base,z,base[2]/base[0],base[3]/base[1],m) for z in (U,W)]
        return dict(count_only=p,count_and_mean_excess=q)
    a,b=[one(f) for f in raw]
    if max(abs(np.asarray(a[k]['conditional_relative_split_interval'])-
               np.asarray(b[k]['conditional_relative_split_interval'])).max() for k in b)>1e-13:
        raise ArithmeticError('Paired error did not resolve the conditional split bounds')
    if abs(b['count_only']['baseline_relative_split']-r['relative_neutrality_split'])>2e-15:
        raise ArithmeticError('Independent stable baseline split disagrees')
    return dict(grid='fine',x=.05,scalar_div=n,charge_div=2,band=r['band'],**b,
                scope='Conservative correlated finite-panel bounds. Fixed pole and charge; moments of each occupied in-band shell constrained. No continuum, physical population, or dynamic error guarantee.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('raw',nargs=2,type=Path)
    p.add_argument('state',type=Path)
    p.add_argument('archive',type=Path)
    p.add_argument('--scalar-div',type=int,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.raw,a.state,a.archive,a.scalar_div),indent=2,allow_nan=False))
