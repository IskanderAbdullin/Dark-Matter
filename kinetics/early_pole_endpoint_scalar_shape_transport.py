"""Measure exact 1-D histogram transport on matched scalar-only endpoint charts."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from early_pole_endpoint_scalar_spectral_moment_budget import load_card


CUT=1e-5


def cdf(d,state,upper):
    e=d.ds.ed
    if not e[0]<upper<=e[-1]:
        raise ArithmeticError('Invalid support')
    k=np.searchsorted(e,upper,side='left')
    edge=np.r_[e[:k],upper]
    f=state[:k]
    l,h=edge[:-1],edge[1:]
    number=f*(h-l)*(h*h+h*l+l*l)
    moment=.75*f*(h-l)*(h**3+h*h*l+h*l*l+l**3)
    total=fsum(number)
    p=number/total
    return edge,np.r_[0.,np.cumsum(p)],np.r_[np.cumsum(p[::-1])[::-1],0.],fsum(moment)/total


def interval(a,b,fa,fb):
    h=b-a
    mean=(1.5*a*a+a*h+.25*h*h)/(3*a*a+3*a*h+h*h)
    return h*(fa+(fb-fa)*mean)


def pair(L,R):
    l,fl,sl,lmu=L
    r,fr,sr,rmu=R
    edge=np.union1d(l,r)
    q3=edge**3
    a,b=np.interp(q3,l**3,fl),np.interp(q3,r**3,fr)
    # In the warm tail, subtract accurately summed upper-tail masses
    # instead of two CDFs already close to one.
    rev=np.interp(q3,l**3,sl)-np.interp(q3,r**3,sr)
    d=np.where((a>.5)&(b>.5),rev,b-a)
    d[0]=d[-1]=0.
    I=[]
    J=[]
    Y=[]
    theta,wt=roots_legendre(8)
    for a,b,z,w in zip(edge[:-1],edge[1:],d[:-1],d[1:],strict=True):
        cuts=[(a,b,z,w)]
        if z*w<0:
            f=-z/(w-z)
            mid=np.cbrt(a**3+f*(b**3-a**3))
            mid=np.clip(mid,a,b)
            cuts=[(a,mid,z,0.),(mid,b,0.,w)]
        for h,k,u,v in cuts:
            if not h<k:
                continue
            val=interval(h,k,u,v)
            Y.append(val)
            I.append(abs(val))
            q=h+(k-h)*(theta+1)/2
            s=(q-h)*(q*q+q*h+h*h)/((k-h)*(k*k+k*h+h*h))
            J.append((k-h)/2*float(wt@abs(u+(v-u)*s)))
    transport=fsum(I)
    signed=fsum(Y)
    if abs(fsum(J)-transport)>4e-13*transport or abs(signed+rmu-lmu)>2e-14*max(lmu,rmu):
        raise ArithmeticError('Histogram transport, exact signed mean, or independent quadrature did not close')
    return dict(wasserstein1_q=transport,total_variation=fsum(abs(np.diff(d)))/2,
                maximum_absolute_cdf_difference=float(max(abs(d))),means=[lmu,rmu],
                signed_mean_shift=rmu-lmu,transport_over_old_mean=transport/lmu,
                absolute_mean_shift_over_transport=abs(rmu-lmu)/transport if transport else None)


def run(paths):
    lev={2:(2,4),3:(1,2,4),4:(1,2,4,8)}[len(paths)]
    ds=[load_card(p,n) for n,p in zip(lev,paths,strict=True)]
    v=[(cdf(d,a,CUT),cdf(d,a,d.ds.ed[-1])) for d,a,_ in ds]
    out=[]
    for i in range(len(lev)-1):
        a,b=v[i:i+2]
        out.append(dict(scalar_levels=[lev[i],lev[i+1]],cold_normalized=pair(a[0],b[0]),full=pair(a[1],b[1])))
    return dict(cold_upper_q=CUT,comparisons=out,scope='One-dimensional transport between numerical histograms; no physical transfer-function inference.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('states',nargs='+',type=Path)
    a=p.parse_args()
    if not 2<=len(a.states)<=4:
        p.error('need S2,S4 or native,S2,S4[,S8]')
    print(json.dumps(run(a.states),indent=2,allow_nan=False))
