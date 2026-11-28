"""Conditional common-phase streaming bounds and constructive atomic witnesses."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_frozen_ballistic_shape import loss
from early_pole_endpoint_radiation_clock_streaming_probe import average_speed
from early_pole_endpoint_scalar_gauss2_energy_reconstruction import shell_rule
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack
from early_pole_endpoint_warm_three_chart_stress_hull import conditioned


def kernel(q,m,A,b):
    return loss(b*average_speed(q,m,A))


def curvature(a,b,m,A,phase):
    low,high=np.hypot([a,b],m)
    h=(b-a)*(b+a)/(high+low)
    # Energy coordinate derivatives for q(xi), exact upper bounds over the shell.
    d1=h*low/a
    d2=h*h*m*m/a**3
    # Uniform in the prescribed clock because its instantaneous mass is >= m.
    cap=(m*A)**2
    v1=min(1./m,1./a,cap/a**3)
    v2=min(3./a**2,3.*b/m**3,3.*cap/a**4)
    # From 1-j0(z)=int_0^1 [1-cos(z*t)] dt.
    return phase**2*(v1*d1)**2/3.+abs(phase)*(v2*d1*d1+v1*d2)/2.


def hull(a,b,m,A,phase,mean,n):
    grid=np.linspace(0.,1.,n+1)
    E=np.hypot(a,m)
    h=(b-a)*(b+a)/(np.hypot(b,m)+E)
    vals=kernel(np.sqrt(a*a+grid*h*(2.*E+grid*h)),m,A,phase)
    i,j=np.flatnonzero(grid<mean),np.flatnonzero(grid>mean)
    zz=[vals[grid==mean]]
    if i.size and j.size:
        i=i[:,None]; j=j[None,:]
        w=(mean-grid[i])/(grid[j]-grid[i])
        zz.append(((1-w)*vals[i]+w*vals[j]).ravel())
    mix=np.concatenate(zz)
    err=curvature(a,b,m,A,phase)/(8*n*n)
    left,right=float(np.min(mix)),float(np.max(mix))
    return np.array([max(0.,left-err),right+err]),np.array([left,right])


def average(chart,m,A,phase,n):
    lo=[]; hi=[]; inner=[]; nodes=[]; w=[]
    for a,b,z,mu,q,p,warm in chart:
        if not z:
            continue
        if warm:
            outer,near=hull(a,b,m,A,phase,mu,n)
            witness=fsum(p*kernel(q,m,A,phase))/z
        else:
            Q,_,L=shell_rule(a,b,m,n)
            W=np.exp(L)
            witness=fsum(W*kernel(Q,m,A,phase))/fsum(W)
            outer=near=np.array([witness,witness])
        if not outer[0]-5e-14<=witness<=outer[1]+5e-14:
            raise ArithmeticError('Nodal witness violates its shell envelope')
        lo.append(z*outer[0]); hi.append(z*outer[1]); inner.append(z*near); nodes.append(z*witness); w.append(z)
    N=fsum(w)
    return np.array([fsum(lo),fsum(hi)])/N,np.array([fsum(v) for v in np.array(inner).T])/N,fsum(nodes)/N


def run(native,accepted,control,stream,state):
    (u,uh)=read(native)
    ab,bb=accepted.read_bytes(),stream.read_bytes()
    a,b=json.loads(ab),json.loads(bb)
    cr=hashlib.sha256(control.read_bytes()).hexdigest()
    ar=hashlib.sha256(ab).hexdigest()
    R=a['rows'][-1]
    C=[read(p) for p in state]
    if (a['charts']!=['G8','G9','mixed300'] or R['x']!=.05 or [h for _,h in C]!=R['states_sha256']
            or b['input_sha256']!=[uh,ar,cr,*R['states_sha256']]):
        raise ArithmeticError('Require exact accepted mature states and prior stream lock')
    d,y=unpack(u,False,.05)
    m=d.state(.05,y)[0]
    if m!=R['common_native_mass'] or b['formal_initial_mass']!=m:
        raise ArithmeticError('Frozen mass cannot be changed')
    C=[conditioned(d,x,m,.05,z) for (x,_),z in zip(C,(None,None,control),strict=True)]
    if [c[0] for c in C]!=[8,9,'mixed300']:
        raise ArithmeticError('Wrong fixed shape charts')
    out=[]
    for O in b['rows']:
        A=O['A']
        q=next(z for z in O['rows'] if z['u']==1.)['representations'][0]
        phase=q['phase'][0]
        coarse,fine=[[average(z,m,A,phase,n) for _,z,_ in C] for n in (128,384)]
        cases=[]
        for i,(g,h,stored) in enumerate(fine):
            if abs(stored-q['matrix'][i][0])>2e-14:
                raise ArithmeticError('Same-phase archived nodal result did not replay')
            if abs(stored-coarse[i][2])>2e-14:
                raise ArithmeticError('Cold constant quadrature did not close')
            c0,c1=coarse[i][0]
            if c0-2e-14>h[0] or c1+2e-14<h[1]:
                raise ArithmeticError('Fine positive witnesses violate the coarser valid envelope')
            cases.append(dict(outer=g.tolist(),constructive=h.tolist(),nodal=stored))
        pairs=[]
        for i,j in ((0,1),(0,2),(1,2)):
            l,u,wi=fine[i]; v,w,wj=fine[j]
            outer=[v[0]-l[1],v[1]-l[0]]
            inner=[w[0]-u[1],w[1]-u[0]]
            if not outer[0]-2e-14<=inner[0]<=inner[1]<=outer[1]+2e-14:
                raise ArithmeticError('Constructive pair bounds should lie inside outer bounds')
            code=('positive' if outer[0]>0 else 'negative' if outer[1]<0 else
                  'both_witnessed' if inner[0]<0<inner[1] else 'inconclusive')
            pairs.append(dict(comparison=[a['charts'][j],a['charts'][i]],difference_outer=outer,
                              difference_constructive=inner,nodal_difference=wj-wi,sign=code))
        out.append(dict(A=A,phase=phase,comparisons=pairs,charts=cases))
    return dict(x=.05,formal_mass=m,rows=out,stream_sha256=hashlib.sha256(bb).hexdigest(),
                scope='Conditional warm positive-measure envelopes at fixed number and mean energy; cold cells held constant and integrated numerically. Outer interpolation remainders evaluated in ordinary floating point; no cosmological transfer claim.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('native','accepted','control','stream'):
        p.add_argument(k,type=Path)
    p.add_argument('state',type=Path,nargs=3)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
