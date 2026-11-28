"""Bound normalized common-q bin profiles under independent finite scalar reweights."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from early_pole_endpoint_cold_panel_reweight_envelope import mean_vertices,panel_excess
from early_pole_endpoint_cold_scalar_band_ladder import verify_pair
from early_pole_endpoint_continuum_cold_band_receipt import read as panels
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart
from early_pole_endpoint_stencil_remap_pilot import SPECS,load


def shares(offset,vertices,base):
    const=fsum(offset)
    if np.min(base)<0 or not const>0 or any(np.min(a)<0 for a in vertices):
        raise ArithmeticError('Need positive finite panel totals')
    totals=[a.sum(axis=1) for a in vertices]
    S=fsum(base)
    results=[]
    for j,b in enumerate(base):
        limits=[]
        for f in (min,max):
            def F(r):
                return (offset[j]-r*const+fsum(f(a[:,j]-r*t)
                         for a,t in zip(vertices,totals,strict=True)))/S
            x=brentq(F,-1e-14,1+1e-14,xtol=4e-16,rtol=1e-15)
            if abs(F(x))>1e-13 or x<-5e-15 or x>1+5e-15:
                raise ArithmeticError('Fractional profile extremum did not resolve')
            limits.append(float(np.clip(x,0,1)))
        if not limits[0]-3e-14<=b/S<=limits[1]+3e-14:
            raise ArithmeticError('The actual profile escaped its allowed reweighting interval')
        results.append(limits)
    p=np.asarray(results)
    if fsum(p[:,0])>1+4e-14 or fsum(p[:,1])<1-4e-14:
        raise ArithmeticError('Invalid probability interval family')
    return p


def fit_state(path,archive,n,pair):
    x,sha=read(path)
    r=verify_pair(pair,path,archive,n)['diagnostic']
    d=ScalarSubshellChart(load('fine')[0],SPECS['fine'],2,n,1e-5)
    return d,np.asarray(x['endpoint_occupation'])[:d.ns],r['full_pole_and_snapshot']['scalar_mass'],sha


def geometry(path):
    u=panels(path)
    e=np.array([u[0]['bounds'][0]]+[x['bounds'][1] for x in u])
    H=np.array([x['coeff'] for x in u])
    return e,H


def build(book,pop,common):
    (e,z),(d,y,m,_)=book,pop
    a,b=e[:-1],e[1:]
    ix=np.searchsorted(d.ds.ed,(a+b)/2)-1
    kk=np.searchsorted(common,(a+b)/2)-1
    if (np.any(ix<0) or np.any(ix>=d.ns) or np.any(kk<0) or np.any(kk>=len(common)-1)
        or np.any((a<d.ds.ed[ix])|(b>d.ds.ed[ix+1]))
        or np.any((a<common[kk])|(b>common[kk+1]))):
        raise ArithmeticError('A q panel crosses an original shell or a shared group boundary')
    f=y[ix]
    V=(b-a)*(b*b+b*a+a*a)/3
    X=panel_excess(a,b,m,48)
    if np.max(abs(X/panel_excess(a,b,m,24)-1))>5e-14:
        raise ArithmeticError('The within-panel conditioning means did not close')
    U=[]; W=[]
    total=len(common)-1
    # Preserve the complete gross number profile, including spontaneous A.
    A=z[:,0,0]/(1+f)
    C=np.bincount(kk,weights=A,minlength=total)
    T=np.bincount(kk,weights=z[:,0,0],minlength=total)
    P=np.zeros((len(f),total))
    P[np.arange(len(f)),kk]=A/V
    if np.max(abs((C+(V*f)@P-T)/np.maximum(T,1e-100)))>4e-14:
        raise ArithmeticError('Original profile fails panel normalization')
    for j in np.unique(ix):
        u=ix==j
        w=V[u]*f[u]
        N=fsum(w)
        if N==0:
            continue
        if len(w)==1:
            mat=np.ones((1,1))
        else:
            mat=mean_vertices(X[u],fsum(w*X[u])/N)
        U.append(N*P[u])
        W.append(N*(mat@P[u]))
    return dict(base=T/fsum(T),count_only=shares(C,U,T),count_and_mean_excess=shares(C,W,T))


def compare(a,b,edges):
    out={}
    for key in ('count_only','count_and_mean_excess'):
        u,v=[z[key] for z in (a,b)]
        gain=np.maximum(0.,v[:,0]-u[:,1])
        loss=np.maximum(0.,u[:,0]-v[:,1])
        j=np.argmax(np.maximum(gain,loss))
        out[key]=dict(robust_total_variation_lower_bound=max(fsum(gain),fsum(loss)),
                      guaranteed_gain_by_common_bin=gain.tolist(),guaranteed_loss_by_common_bin=loss.tolist(),
                      biggest_certified_difference_q=edges[j:j+2].tolist(),
                      left_right_allowed_share_intervals=[u.tolist(),v.tolist()])
    return out


def run(a,b,states,archive,l):
    if tuple(l) not in ((2,4),(4,8),(8,16)):
        raise ArithmeticError('Use consecutive accepted scalar levels')
    qa,qb=[fit_state(p,archive,n,files) for p,n,files in zip(states,l,(a,b),strict=True)]
    outs=[]
    for j in (0,1):
        A,B=map(geometry,(a[j],b[j]))
        e=np.intersect1d(A[0],B[0])
        if len(e)<5 or [e[0],e[-1]]!=[1e-6,1.5e-5]:
            raise ArithmeticError('Need the complete shared cold support')
        u,v=build(A,qa,e),build(B,qb,e)
        outs.append((e.tolist(),[u['base'].tolist(),v['base'].tolist()],compare(u,v,e)))
    (e,a,u),(f,b,v)=outs
    if e!=f or np.max(abs(np.array(a)-np.array(b)))>1e-8:
        raise ArithmeticError('Both panel meshes must agree')
    err=max(np.max(abs(np.asarray(u[k]['left_right_allowed_share_intervals'])-
               np.asarray(v[k]['left_right_allowed_share_intervals']))) for k in v)
    if err>1e-8:
        raise ArithmeticError('Coarse and fine profile interval witnesses do not agree')
    return dict(x=.05,grid='fine',charge_div=2,scalar_levels=l,state_sha256_uncompressed=[qa[-1],qb[-1]],
                shared_q_edges=e,base_normalized_profiles=b,**v,maximum_pair_mesh_absolute_interval_difference=err,
                scope='Finite positive panel reweighting at each stored chart own optical fit, whole band normalization, no inference within a panel or beyond the computed supports.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--left',nargs=2,type=Path,required=True)
    p.add_argument('--right',nargs=2,type=Path,required=True)
    p.add_argument('--states',nargs=2,type=Path,required=True)
    p.add_argument('--levels',nargs=2,type=int,required=True)
    p.add_argument('--archive',type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.left,a.right,a.states,a.archive,a.levels),indent=2,allow_nan=False))
