"""Finite-panel envelope at fixed scalar cell counts, optionally also mean excess."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import roots_legendre

from early_pole_endpoint_cold_scalar_band_ladder import verify_pair
from early_pole_endpoint_continuum_cold_band_receipt import read as raw_panels
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart
from early_pole_endpoint_stencil_remap_pilot import SPECS,load


def panel_excess(a,b,m,n):
    t,w=roots_legendre(n)
    q=a[:,None]+(b-a)[:,None]*(t+1)/2
    z=w*q*q
    return np.sum(z*q*q/(np.hypot(q,m)+m),axis=1)/np.sum(z,axis=1)


def mean_vertices(x,target):
    """Weights at vertices of a unit simplex with one fixed mean."""
    one=[]
    for i,u in enumerate(x):
        if u==target:
            z=np.zeros(len(x))
            z[i]=1.
            one.append(z)
        if u>=target:
            continue
        for j in np.flatnonzero(x>target):
            z=np.zeros(len(x))
            t=(target-u)/(x[j]-u)
            z[[i,j]]=[1-t,t]
            one.append(z)
    a=np.asarray(one)
    if (not len(a) or np.min(a)<0 or np.max(abs(a.sum(axis=1)-1))>2e-15
            or np.max(abs((a@x)/target-1))>2e-14):
        raise ArithmeticError('Conditional positive panel vertices did not close')
    return a


def envelope(const,arr,base):
    # Every row is the total stirred contribution at one cell vertex.
    # For each independent shell, summing independent extrema is exact.
    P=[const+np.array([fsum(f(x[:,j]) for x in arr) for j in range(3)]) for f in (min,max)]
    if min(P[0])<=0 or np.any(P[0]>base) and np.max(P[0]/base-1)>2e-14 or np.max(base/P[1]-1)>2e-14:
        raise ArithmeticError('Baseline coefficient left its finite-panel interval')
    # Correlate B and D per chosen redistribution, rather than dividing
    # independently optimized envelope components.
    roots=[]
    anchor=base[1]/base[2]
    bracket=[0.,2*float(P[1][1]/P[0][2]/anchor)+1.]
    for fn in (min,max):
        def z(r):
            n=r*anchor
            return (const[1]-n*const[2]+fsum(fn(x[:,1]-n*x[:,2]) for x in arr))/base[1]
        y=brentq(z,*bracket,xtol=2e-14)
        if abs(z(y))>1e-12:
            raise ArithmeticError('Conditional fractional envelope not resolved')
        roots.append(y)
    if not roots[0]-1e-13<=1<=roots[1]+1e-13:
        raise ArithmeticError('Baseline neutral root left its independent panel hull')
    return dict(coefficient_relative_to_original_min_max=[(x/base).tolist() for x in P],
                conditional_neutrality_relative_to_original_min_max=roots)


def one(path,d,y,m):
    raw=raw_panels(path)
    e=np.array([raw[0]['bounds'][0]]+[z['bounds'][1] for z in raw])
    a,b=e[:-1],e[1:]
    ir=np.searchsorted(d.ds.ed,(a+b)/2)-1
    if np.any(ir<0) or np.any(ir>=d.ns) or np.any(a<d.ds.ed[ir]) or np.any(b>d.ds.ed[ir+1]):
        raise ArithmeticError('A q panel crosses its archive scalar cell')
    f=y[ir]
    V=(b-a)*(b*b+b*a+a*a)/3
    X=panel_excess(a,b,m,48)
    if np.max(abs(X/panel_excess(a,b,m,24)-1))>5e-14:
        raise ArithmeticError('Within-panel energy means did not resolve')
    z=np.array([u['coeff'][0] for u in raw])
    A=z[:,0]/(1+f)
    B=np.divide(z[:,1],f,out=np.zeros_like(f),where=f>0)
    if np.any(B>A) or np.any((f==0)&(z[:,1]!=0)):
        raise ArithmeticError('This finite positive envelope requires a positive affine slope')
    coeff=np.stack((A,B,A-B),axis=-1)/V[:,None]
    const=np.array([fsum(A),0.,fsum(A)])
    base=np.array([fsum(v) for v in z.T])
    if np.max(abs((const+(f*V)@coeff)/base-1))>5e-14:
        raise ArithmeticError('Original finite scalar step source did not close')
    summary=[]
    p0=[]; p1=[]
    for j in np.unique(ir):
        ix=ir==j
        p=V[ix]*f[ix]
        u=fsum(p)
        s=coeff[ix]
        if u==0:
            continue
        t=fsum(p*X[ix])/u
        if len(p)==1:
            opt=s
        else:
            opt=mean_vertices(X[ix],t)@s
        p0.append(u*s)
        p1.append(u*opt)
        summary.append((int(j),int(len(p)),int(len(opt))))
    return dict(count_only=envelope(const,p0,base),count_and_mean_excess=envelope(const,p1,base),
                conditioned_nonempty_cells=len(summary),max_panels_in_a_cell=max(z[1] for z in summary),
                max_two_constraint_vertices_in_a_cell=max(z[2] for z in summary))


def run(paths,state,archive,n):
    if n not in (2,4,8,16):
        raise ArithmeticError('Need an independently resolved cold-only scalar refinement')
    r=verify_pair(paths,state,archive,n)['diagnostic']
    a,_=read(state)
    y=np.asarray(a['endpoint_occupation'])
    d=ScalarSubshellChart(load('fine')[0],SPECS['fine'],2,n,1e-5)
    m=r['full_pole_and_snapshot']['scalar_mass']
    lo,hi=[one(f,d,y,m) for f in paths]
    err=max(np.max(abs(np.asarray(lo[k][j])-np.asarray(hi[k][j])))
            for k in ('count_only','count_and_mean_excess')
            for j in hi[k])
    if err>4e-8:
        raise ArithmeticError('Finite reweighting bounds changed between accepted panels')
    return dict(grid='fine',x=.05,charge_div=2,scalar_div=n,band=r['band'],columns=('forward','inverse','slope'),
                fixed_scalar_mass=m,**hi,maximum_absolute_paired_mesh_defect=err,
                scope='Finite piecewise-constant q-panel reweights on retained support, fixing each shell amount wholly inside the band and optionally its mean E-m. Charges and pole held. Not a continuum extremum, actual evolving uncertainty or photon prediction.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('raw',nargs=2,type=Path)
    p.add_argument('state',type=Path)
    p.add_argument('archive',type=Path)
    p.add_argument('--scalar-div',type=int,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.raw,a.state,a.archive,a.scalar_div),indent=2,allow_nan=False))
