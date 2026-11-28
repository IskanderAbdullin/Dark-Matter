"""Batched fixed-aperture spinors with exact two-arm local histogram support."""
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import roots_legendre

from early_pole_hard_yukawa_photoproduction_reference import I, dot, norm_spinors, slash
from finite_mass_qed_frozen_transverse_plasmon_gate_reference import Q, Species
from finite_mass_qed_longitudinal_static_spectral_reference import DensityResponse
from finite_mass_qed_timelike_pole_affine_source_reference import run_one


def pole(ch,name,k):
    if name=='L':
        return DensityResponse(ch).pole(k)
    def eq(w):
        return (w-k)*(w+k)-sum(s.field(k/w)[0] for s in ch)
    w=brentq(eq,np.nextafter(k,np.inf),k+.006,xtol=5e-17)
    return w,1/(2*w+k/w**2*sum(s.field(k/w)[1] for s in ch))


def reduce_vertex(p,b,K,m0,m1,e):
    # All momenta are contravariant; spinors keep the full angular layout.
    w,k=K[0],K[3]
    mu2=(w-k)*(w+k)
    v=[]
    for o,m in ((p,m0),(b,m1)):
        spatial=o[...,1:]
        z=spatial[...,2]
        q=np.linalg.norm(spatial,axis=-1)
        turn=q-z
        mask=z>0
        turn[mask]=np.sum(spatial[mask,:2]**2,axis=-1)/(q[mask]+z[mask])
        v.append((m*m*w*w+mu2*q*q)/(np.hypot(m,q)*w+q*k)+k*turn)
    left=(2*dot(p,e)[...,None,None]*I+slash(K)@slash(e))/(mu2+2*v[0])[...,None,None]
    right=(2*dot(b,e)[...,None,None]*I-slash(e)@slash(K))/(mu2-2*v[1])[...,None,None]
    return left+right


def full_vertex(p,b,K,m0,m1,e):
    s=dot(K,K)
    return ((slash(p+K)+m0*I)@slash(e)/(s+2*dot(p,K))[...,None,None]
             +slash(e)@(slash(b-K)+m1*I)/(s-2*dot(b,K))[...,None,None])


def radial(ch,k,w,R,es,c,name,n,az,audit=False):
    P,F=ch
    M,m=P.m,F.m
    T=np.array([es-w,Q*np.sqrt((1-c)*(1+c)),0.,Q*c-k])
    v=np.linalg.norm(T[1:])
    t=T[0]
    s=(t-v)*(t+v)
    if t>=0:
        raise ValueError('This branch requires a past-directed transfer energy')
    if s>=0:
        return np.zeros(3),False,0.
    A=M*M+s-m*m
    lo=(-A*t+v*np.sqrt(A*A-4*M*M*s))/(-2*s)
    hi=max(P.E[-1],F.E[-1]+t)
    if not lo<hi:
        return np.zeros(3),False,0.
    ee=np.unique(np.clip(np.r_[lo,hi,P.E,F.E+t],lo,hi))
    a,b=ee[:-1],ee[1:]
    c0=(a+b)/2
    ff=np.r_[0.,P.f,0.][np.searchsorted(P.E,c0)]
    gg=np.r_[0.,F.f,0.][np.searchsorted(F.E,c0-t)]
    keep=(ff>0)|(gg>0)
    if not any(keep):
        return np.zeros(3),False,0.
    a,b,ff,gg,c0=(x[keep] for x in (a,b,ff,gg,c0))
    x,h=roots_legendre(n)
    e=(c0[:,None]+(b-a)[:,None]*x/2).ravel()
    z=(e*t-A/2)/v
    p2=(e-M)*(e+M)
    radicand=p2-z*z
    if np.min(radicand/(1+p2))<-3e-14:
        raise ArithmeticError('No valid radial angle')
    p=np.sqrt(np.maximum(radicand,0.))
    axis=T[1:]/v
    cross=np.array([axis[2],0.,-axis[0]])
    aphi=(np.arange(az)+.5)*2*np.pi/az
    transverse=np.cos(aphi)[:,None]*cross+np.sin(aphi)[:,None]*np.array([0.,1.,0.])
    mom=z[:,None,None]*axis+p[:,None,None]*transverse[None,:,:]
    P0=np.concatenate([np.broadcast_to(e[:,None,None],(*mom.shape[:2],1)),mom],axis=-1)
    B0=P0-T
    K=np.array([w,0.,0.,k])
    pol=np.eye(4)[1:3] if name=='T' else np.eye(4)[:1]
    integral=0.
    max_def=0.
    for u in pol:
        A=norm_spinors(P0,B0,M,m,reduce_vertex(P0,B0,K,M,m,u))
        if audit:
            B=norm_spinors(P0,B0,M,m,full_vertex(P0,B0,K,M,m,u))
            max_def=max(max_def,float(np.sum(abs(A-B))/np.sum(B)))
        integral=integral+A.mean(axis=-1)
    e_int=((b-a)/2)*((integral.reshape(-1,n))@h)*R/v
    C=np.asarray([ff*(1-gg),gg*(1-ff),ff-gg])@e_int
    if min(C[:2])<0 or abs(C[0]-C[1]-C[2])>5e-13*max(C[:2]):
        raise ArithmeticError('Impossible signed population support')
    return C,True,max_def


def load(g):
    folder=Path(__file__).parent/'results'
    with gzip.open(folder/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as s:d=json.load(s)
    M=d['x_end']
    ch=[Species(m,d[f'endpoint_{j}_momentum_edges'],d[f'endpoint_{j}_occupations'])
        for m,j in ((M,'parent'),(M*d['mass_ratio'],'fermion'))]
    es=np.hypot(Q,d['checkpoints'][-1]['scalar_pole_over_T'])
    return d,ch,es


def fixture():
    out=[]
    for g in ('fine','third'):
        d,ch,es=load(g)
        for name,k in (('T',.0015),('T',.00256),('L',.0012),('L',.00256)):
            w,R=pole(ch,name,k)
            for c in (0.,-.77,.82):
                A,ok,trace=radial(ch,k,w,R,es,c,name,3,6,True)
                B=run_one(d,ch,g,name,k,c,3,6)['raw_kernel']
                err=float(max(abs(A/B-1.)))
                if not ok or max(err,trace)>2e-7:
                    raise ArithmeticError((g,name,k,c,err,trace))
                out.append(dict(grid=g,pol=name,k=k,c=c,batch_vs_fixed_rel=err,vertex_trace_rel=trace))
    return out


if __name__=='__main__':
    print(json.dumps(fixture(),indent=2,allow_nan=False))
