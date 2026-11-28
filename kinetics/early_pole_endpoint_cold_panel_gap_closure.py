"""Finite-panel self-consistent frozen Hartree shifts; no re-evaluated collision rates."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import roots_legendre

from early_pole_endpoint_cold_panel_reweight_envelope import mean_vertices,panel_excess
from early_pole_endpoint_cold_scalar_band_ladder import verify_pair
from early_pole_endpoint_continuum_cold_band_receipt import read as raw_panels
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart
from early_pole_endpoint_stencil_remap_pilot import SPECS,load


def nodes(a,b,n):
    a,b=np.asarray(a),np.asarray(b)
    x,w=roots_legendre(n)
    q=a[:,None]+(b-a)[:,None]*(x+1)/2
    v=(b-a)[:,None]*q*q*w/(4*np.pi**2)
    V=(b-a)*(b*b+b*a+a*a)/(6*np.pi**2)
    if np.max(abs(v.sum(axis=1)/V-1))>3e-14:
        raise ArithmeticError('Positive volume moment mismatch')
    return q,v,V


def one(d,f,m,edges,n):
    x=.05
    start,end=edges[0],edges[-1]
    a,b=edges[:-1],edges[1:]
    q,v,V=nodes(a,b,n)
    outside=np.union1d(d.ds.ed,[start,end])
    A,B=outside[:-1],outside[1:]
    stay=(B<=start)|(A>=end)
    A,B=A[stay],B[stay]
    z,w,Vo=nodes(A,B,n)
    keep=np.searchsorted(d.ds.ed,(A+B)/2)-1
    ir=np.searchsorted(d.ds.ed,(a+b)/2)-1
    if np.any((ir<0)|(ir>=d.ns)) or np.any((a<d.ds.ed[ir])|(b>d.ds.ed[ir+1])):
        raise ArithmeticError('A finite panel crossed a declared scalar shell')
    ref=V*f[ir]
    fixed=Vo*f[keep]
    if abs((fsum(ref)+fsum(fixed))/fsum(d.ds.w*f)-1)>3e-14:
        raise ArithmeticError('Cold and fixed parts do not preserve the scalar count')
    oz=(w*f[keep,None]).ravel()
    ozq=z.ravel()
    m0=1e-9*x
    lam=d.ds.lam
    if lam<=0:
        raise ArithmeticError('Need the positive Hartree coupling')
    nodal=fsum(d.ds.w*f/np.hypot(d.ds.p,m))
    if abs((m*m-m0*m0)/(lam/2*nodal)-1)>2e-11:
        raise ArithmeticError('Reported pole is not the nodal gap root')
    h=panel_excess(a,b,m,n)
    qids=[]
    for j in np.unique(ir):
        i=ir==j
        N=fsum(ref[i])
        if N<=0:
            continue
        E=fsum(ref[i]*h[i])/N
        weights=mean_vertices(h[i],E) if np.count_nonzero(i)>1 else np.ones((1,1))
        qids.append((i,N,weights))
    def row(M,mode,fun):
        outside=fsum(oz/np.hypot(ozq,M))
        f0=np.sum(v/np.hypot(q,M),axis=1)/V
        inside=fsum(ref*f0) if mode==0 else fsum(
            N*fun(f0[i] if mode==1 else (w@f0[i])) for i,N,w in qids)
        return (M-m0)*(M+m0)-lam/2*(outside+inside)
    lim=[]
    for mode in (0,1,2):
        r=[]
        for fun in (min,max):
            a=brentq(lambda z:row(z,mode,fun),.8*m,1.2*m,xtol=1e-24,rtol=1e-15)
            if abs(row(a,mode,fun))>1e-13*m*m:
                raise ArithmeticError('The finite Hartree root solve did not close')
            r.append(float(a))
        lim.append(r)
    a,c,e=[np.array(z) for z in lim]
    if abs(a[0]/a[1]-1)>1e-15 or not c[0]*(1-3e-15)<=e[0]<=a[0]*(1+3e-15) or not a[0]*(1-3e-15)<=e[1]<=c[1]*(1+3e-15):
        raise ArithmeticError('The baseline or narrower conditional root escaped the hull')
    return dict(base_step_continuum_root_relative_to_original_nodal_minus_one=a[0]/m-1,
                count_only_root_min_max_relative_to_nodal_minus_one=(c/m-1).tolist(),
                count_and_mean_excess_root_min_max_relative_to_nodal_minus_one=(e/m-1).tolist(),
                count_only_root_width_relative_to_nodal=(c[1]-c[0])/m,
                count_and_mean_excess_root_width_relative_to_nodal=(e[1]-e[0])/m)


def run(raw,state,archive,n):
    r=verify_pair(raw,state,archive,n)
    p=raw_panels(raw[1])
    edges=np.array([p[0]['bounds'][0]]+[v['bounds'][1] for v in p])
    card,_=read(state)
    d=ScalarSubshellChart(load('fine')[0],SPECS['fine'],2,n,1e-5)
    m=r['diagnostic']['full_pole_and_snapshot']['scalar_mass']
    f=np.asarray(card['endpoint_occupation'])[:d.ns]
    lo,hi=[one(d,f,m,edges,n) for n in (24,48)]
    err=max(np.max(abs(np.asarray(lo[k])-np.asarray(hi[k]))) for k in hi)
    if err>1e-13:
        raise ArithmeticError('Finite-envelope root comparison did not resolve')
    return dict(grid='fine',scalar_div=n,x=.05,charge_div=2,nodal_gap=m,band=edges[[0,-1]].tolist(),**hi,
                numerical_absolute_difference_between_meshes=err,
                scope='Finite panel reweighting within the retained band, fixed original out-of-band shape and frozen Hartree closure. No re-evaluated cold rates, true continuum extremum or dynamic physical uncertainty.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('raw',type=Path,nargs=2)
    p.add_argument('state',type=Path)
    p.add_argument('archive',type=Path)
    p.add_argument('--scalar-div',type=int,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.raw,a.state,a.archive,a.scalar_div),indent=2,allow_nan=False))
