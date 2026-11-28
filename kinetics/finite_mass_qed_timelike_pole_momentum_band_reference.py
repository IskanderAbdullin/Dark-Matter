"""Bounded timelike pole aperture with exact angular/end-cell momentum joins."""
import argparse
from functools import lru_cache
import json
import os
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import roots_legendre

from finite_mass_qed_frozen_transverse_plasmon_gate_reference import Q
from finite_mass_qed_timelike_pole_angular_reference import angle
from finite_mass_qed_timelike_pole_batched_absorption_reference import load, pole


@lru_cache(None)
def mesh(n):
    return roots_legendre(n)


def boundaries(ch,state,es,limits,samples):
    M,m=[x.m for x in ch]
    def bisect_value(k,E,index,cos):
        t=es-state(k)[0]
        if index:
            z=E+t; y=E
        else:
            z=E; y=E-t
        good=(z>=M)&(y>=m)
        pp=np.sqrt(np.maximum(0.,(z-M)*(z+M)))
        pb=np.sqrt(np.maximum(0.,(y-m)*(y+m)))
        rr=pp+pb
        sep=np.divide(M*M-m*m-2*z*t+t*t,rr,out=np.full_like(z,np.nan,dtype=float),where=good&(rr>0))
        return k-Q*cos-sep
    u=np.linspace(*limits,samples)
    ret=[]
    for sign in (-1.,1.):
        phi=lambda k:k-Q*sign+es-state(k)[0]
        vals=[phi(x) for x in u]
        ret.extend(float(x) for x,v in zip(u,vals,strict=True) if v==0.)
        for a,b,v,w in zip(u[:-1],u[1:],vals[:-1],vals[1:],strict=True):
            if v*w<0:
                ret.append(float(brentq(phi,a,b,xtol=5e-17)))
        for i,c in enumerate(ch):
            rr=np.array([bisect_value(k,c.E,i,sign) for k in u])
            ret.extend(float(u[j]) for j in np.flatnonzero(np.any(rr==0.,axis=1)))
            for cell,edge in zip(*np.nonzero(rr[:-1]*rr[1:]<0),strict=True):
                f=lambda k:float(bisect_value(k,c.E[edge],i,sign))
                ret.append(float(brentq(f,u[cell],u[cell+1],xtol=5e-17)))
    return np.unique(np.r_[limits,ret])


def run(grid,name,refine,path):
    d,ch,es=load(grid)
    state=lru_cache(None)(lambda k:pole(ch,name,k))
    limits=(.0012,.003) if name=='T' else (.00085,.003)
    p=boundaries(ch,state,es,limits,49)
    q=boundaries(ch,state,es,limits,89)
    if p.shape!=q.shape or np.max(abs(p-q))>7e-16:
        raise ArithmeticError(('Momentum topology capture',grid,name,len(p),len(q)))
    # Use the same deterministic boundaries in both replay orders.
    bins=q
    outer,cos,inner=(3,6,(6,12)) if refine else (2,3,(3,8))
    x,w=mesh(outer)
    old={}
    if path.is_file():
        for line in path.read_text().splitlines():
            rec=json.loads(line)
            i=rec['i']
            if i in old or rec['range']!=bins[i:i+2].tolist() or rec['mesh']!=[outer,cos,*inner]:
                raise ArithmeticError('Invalid checkpoint panel')
            old[i]=rec
    with path.open('a',buffering=1) as stream:
        for i,(a,b) in enumerate(zip(bins[:-1],bins[1:],strict=True)):
            if i in old:
                continue
            k=(a+b)/2+(b-a)*x/2
            vals=[]; count=0; err=0.
            for z in k:
                A,C,_,D=angle(ch,z,*state(z),es,name,cos,inner)
                vals.append(A*z*z)
                count+=int(C[2])
                err=max(err,D)
            val=((b-a)/2*(w@np.asarray(vals))).tolist()
            print(json.dumps(dict(grid=grid,pol=name,mesh=[outer,cos,*inner],i=i,total=len(bins)-1,range=[a,b],
                                  raw_signed_coefficients=val,forward_only_angle_nodes=count,
                                  max_exact_angular_knot_defect=err),allow_nan=False),file=stream,flush=True)
            os.fsync(stream.fileno())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=('fine','third'))
    p.add_argument('pol',choices=('T','L'))
    p.add_argument('path',type=Path)
    p.add_argument('--refine',action='store_true')
    o=p.parse_args()
    run(o.grid,o.pol,o.refine,o.path)
