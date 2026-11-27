"""Classical timelike pole complement to a closed-path spacelike HTL test."""
import argparse
from functools import lru_cache
import json

import numpy as np
from scipy.optimize import brentq
from early_pole_finite_wilson_htl_reference import gauss, seg


@lru_cache(None)
def modes(k):
    # Solve in log(x-1) to retain the longitudinal residue near the light cone.
    def pieces(u):
        d=np.exp(u)
        x=1+d
        b=d*(d+2)
        h=.5*np.log1p(2/d)
        return x,b,h
    def eqL(u):
        x,b,h=pieces(u)
        return k*k+1-x*h
    def eqT(u):
        x,b,h=pieces(u)
        return -k*k*b+.5*(x*x-x*b*h)
    if k<4.:
        ul=brentq(eqL,-100.,max(8.,np.log(3/k)))
        l,bl,hl=pieces(ul)
        # d[x atanh(1/x)]/dx = atanh(1/x)-x/(x²-1).
        aL=2/l/(l/bl-hl)
    else:
        l,aL=1.,0.
    ut=brentq(eqT,-100.,max(8.,np.log(3/k)))
    x,b,h=pieces(ut)
    # Stable derivative Pi_T'=1.5 x-(3x²-1)h/2.
    aT=2/x/abs(-2*k*k*x+1.5*x-(3*x*x-1)*h/2)
    assert min(aL,aT)>=0 and x>1 and abs(eqT(ut))<1e-4
    return l,aL,x,aT


def js(k,X,mu,t,r,v,mode):
    if mode=='rectangle':
        sine=seg(k*t*X,k*t*X/2)
        qu=k*mu*r
        J0=t*sine*(-np.expm1(-1j*qu))
        Jz=r*seg(-qu,-qu/2)*np.expm1(1j*k*t*X)
    else:
        a=k*t*(X-v*mu)
        b=a-k*mu*r
        A=seg(a,a/2)
        B=seg(b,b/2)
        J0=t*(A-B)
        Jz=v*J0+r*(seg(-k*mu*r,(a+b)/2)-B)
    return np.array([abs(J0)**2,abs(Jz)**2*(1-mu*mu)])


def pole_fold(t,r,v,mode,nrad,nangle,kmax):
    a,w=gauss(nangle)
    # Both signs of omega represented via parity on the full mu range.
    grid=np.r_[-1.,0.,1-np.geomspace(.001,.5,10)[::-1],1.]
    mu=((grid[1:]+grid[:-1])[:,None]/2+np.diff(grid)[:,None]*a/2).flatten()
    mw=(np.diff(grid)[:,None]*w/2).flatten()
    ar,wr=gauss(nrad)
    g=np.r_[0.,np.geomspace(.002,kmax,70)]
    kr=(g[1:]+g[:-1])[:,None]/2+np.diff(g)[:,None]*ar/2
    kw=(np.diff(g)[:,None]*wr/2).flatten()
    rows=[]
    for k in kr.flatten():
        L,aL,T,aT=modes(k)
        wL=aL*(js(k,L,mu,t,r,v,mode)[0]@mw)
        wT=aT*(js(k,T,mu,t,r,v,mode)[1]@mw)
        w0=(js(k,1.,mu,t,r,v,mode)[1]@mw)/(k*k)
        rows.append(np.array([wL,wT,w0])*k*k/(2*np.pi*t))
    return kw@np.array(rows)


def low_k_check(t=995.734, r=.60102770, v=.99939):
    k=1e-7
    x=1/(np.sqrt(3)*k)
    lam=t/np.sqrt(3)
    f=1-2*np.sin(lam)/lam+2*(1-np.cos(lam))/(lam*lam)
    m,w=gauss(100)
    s=js(k,x,m,t,r,v,'triangle')
    # J0 is O(k) at fixed plasmon frequency, which cancels the
    # small-k enhancement in the longitudinal Coulomb pole.
    longitudinal=.5*(w@s[0])*x*x/t**2
    transverse=.5*(w@s[1])/t**2
    ref=(r/t)**2*f/3
    assert abs(longitudinal/ref-1)<1e-5
    assert abs(transverse/(2*ref)-1)<1e-5
    small_omega=np.array([np.max(js(q,.5,m,t,r,v,'triangle')) for q in (1e-8,2e-8,4e-8)])
    assert max(abs(small_omega[1:]/small_omega[:-1]/4-1))<1e-6
    eps=1e-3
    l,_,s,_=modes(eps)
    assert max(abs(np.array([l,s])*eps-1/np.sqrt(3)))<2e-6
    return dict(k=k,longitudinal_projected=longitudinal,transverse_projected=transverse,
                analytic_longitudinal=ref,spacelike_max_power=small_omega.tolist())


def run(level):
    nr,na=(6,20) if level==0 else (12,38) if level==1 else (22,68)
    rows=[]
    for kt in (8.,20.,40.):
        for mode,t,r,v in [('rectangle',32.,.60103446,0.),('rectangle',160.,.60103446,0.),
                           ('triangle',32.,.60103446,.97),('triangle',160.,.60103446,.994),
                           ('triangle',500.,.60103446,.998),('triangle',996.,.60103446,.9982),
                           ('triangle',32.,0.,.97)]:
            x=pole_fold(t,r,v,mode,nr,na,kt)
            assert min(x)>=0.
            if r==0.: assert max(x)==0.
            rows.append(dict(mode=mode,t=t,final_r=r,v=v,kmax=kt,L_pole=x[0],
                             T_pole=x[1],T_free_classical=x[2],
                             T_difference=x[1]-x[2]))
    return dict(scope='Classical HTL timelike residues on a finite test contour, '
                'same-cutoff free pole comparator, not a microscopic linewidth.',
                level=level,low_k_validation=low_k_check(),rows=rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--level',type=int,default=0)
    print(json.dumps(run(p.parse_args().level),indent=2,allow_nan=False))
