"""Finite cold, gap, and shoulder scalar q-bands at one transverse pole momentum."""
import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import roots_legendre

from finite_mass_qed_timelike_pole_angular_reference import angle
from finite_mass_qed_timelike_pole_batched_absorption_reference import load,pole


def q_knots(d,ch,k,w,ms,limits,n):
    M,m=[x.m for x in ch]
    def eval_(q,c,i,E):
        t=np.hypot(ms,q)-w
        x,y=(E+t,E) if i else (E,E-t)
        p=np.sqrt(np.maximum((x-M)*(x+M),0.))
        b=np.sqrt(np.maximum((y-m)*(y+m),0.))
        val=np.divide(M*M-m*m-2*x*t+t*t,p+b,out=np.full_like(x,np.nan,dtype=float),
                      where=(x>=M)&(y>=m)&((p+b)>0.))
        return k-q*c-val
    grid=np.linspace(*limits,n)
    out=[]
    for c in (-1.,1.):
        f=lambda q:k-q*c+np.hypot(ms,q)-w
        vals=[f(x) for x in grid]
        for a,b,v,z in zip(grid[:-1],grid[1:],vals[:-1],vals[1:],strict=True):
            if v*z<0.:
                out.append(float(brentq(f,a,b,xtol=5e-18)))
        for i,x in enumerate(ch):
            V=np.array([eval_(q,c,i,x.E) for q in grid])
            for j,p in zip(*np.nonzero(V[:-1]*V[1:]<0),strict=True):
                u=lambda q:float(eval_(q,c,i,x.E[p]))
                out.append(float(brentq(u,grid[j],grid[j+1],xtol=5e-18)))
    q=np.asarray(d['scalar_momentum_edges'])
    sel=q[(q>limits[0])&(q<limits[1])]
    return np.unique(np.r_[limits,out,sel])


REGIONS={'cold':(1e-6,1.5e-5),'gap':(1.5e-5,4e-5),'shoulder':(4e-5,3e-4)}


def run(g,m,region,output):
    d,ch,_=load(g)
    k=.0015
    w,R=pole(ch,'T',k)
    ms=d['checkpoints'][-1]['scalar_pole_over_T']
    interval=REGIONS[region]
    a=q_knots(d,ch,k,w,ms,interval,49)
    b=q_knots(d,ch,k,w,ms,interval,89)
    if a.shape!=b.shape or max(abs(a-b))>9e-16:
        raise ArithmeticError(('Unstable scalar energy intersections',g,len(a),len(b)))
    n,c,E,Z=((2,3,3,8),(4,6,6,12),(8,8,8,16))[m]
    p,x=roots_legendre(n)
    Q=np.asarray(d['scalar_momentum_edges'])
    F=np.r_[0.,d['endpoint_scalar_occupations'],0.]
    old={}
    if output.is_file():
        for row in map(json.loads,output.read_text().splitlines()):
            j=row['panel']
            if (j in old or not 0<=j<len(b)-1 or row['grid']!=g or row['total']!=len(b)-1
                    or row['limits']!=b[j:j+2].tolist() or row['region']!=region or row['mesh']!=[n,c,E,Z]):
                raise ArithmeticError('Invalid partial q-band receipt')
            old[j]=row
    with output.open('a',buffering=1) as f:
        for j,(l,h) in enumerate(zip(b[:-1],b[1:],strict=True)):
            if j in old:continue
            data=np.zeros((3,3))
            for q,wt in zip((l+h)/2+(h-l)*p/2,(h-l)*x/2,strict=True):
                Es=np.hypot(ms,q)
                occ=F[np.searchsorted(Q,q)]
                A,B,D=angle(ch,k,w,R,Es,'T',c,(E,Z),q)[0]
                y=np.array([(1+occ)*A,occ*B,occ*D+A])
                v=wt*q*q/Es*y
                data[0]+=v
                data[1]+=Es*v
                data[2]+=(q*q/(Es+ms))*v
            f.write(json.dumps(dict(grid=g,panel=j,total=len(b)-1,limits=[l,h],mesh=[n,c,E,Z],
                                    region=region,number_scalar_energy_and_centered_coefficients=data.tolist()),
                               allow_nan=False)+'\n')
            f.flush()
            os.fsync(f.fileno())


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('grid',choices=('fine','third'))
    p.add_argument('mesh',type=int,choices=(0,1,2))
    p.add_argument('region',choices=tuple(REGIONS))
    p.add_argument('out',type=Path)
    o=p.parse_args()
    run(o.grid,o.mesh,o.region,o.out)
