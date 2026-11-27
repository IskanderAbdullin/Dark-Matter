"""Finite Abelian closed-current spacelike HTL covariance benchmark."""
from functools import lru_cache
import argparse
import json

import numpy as np
from scipy.special import roots_legendre

from early_pole_static_longitudinal_dipole_reference import static_shape
from early_pole_moving_soft_dipole_reference import weights


@lru_cache(None)
def gauss(n):
    return roots_legendre(n)


def seg(angle,mid):
    return np.sinc(angle/(2*np.pi))*np.exp(1j*mid)


def radial_loop(k, mode, t, r, v, nangle, nfreq, rx=0., az=0.):
    nodes, aw = gauss(nangle)
    mu = (nodes+1)/2
    fm, fw = gauss(nfreq)
    mu = mu[:, None, None]
    base = np.array([-1.,-.99,-.9,-.6,-.25,0,.25,.6,.9,.99,1.])
    if mode=='rectangle':
        assert rx==0.
    dr = np.sqrt(1-mu*mu)*np.cos(az)
    centers = [0.*mu] if mode=='rectangle' else [v*mu,(v+r/t)*mu+(rx/t)*dr]
    spans = [np.broadcast_to(base[:,None],(base.size,nangle))]
    # Spectral edge and the nearly magnetostatic cut for low k.
    scale=k*k
    spans.append(np.array([-.1,-.01,-.001,-.0001,0.,.0001,.001,.01,.1])[:,None]
                 *np.ones((1,nangle)))
    spans.append(np.array([-100.,-10.,-1.,-.1,.1,1.,10.,100.])[:,None]
                 *scale*np.ones((1,nangle)))
    for c in centers:
        # Localize the two finite-time peaks and logarithmic sinc wings.
        spans.append(c.reshape(-1) + np.array([-64.,-16.,-4.,-1.,0.,1.,4.,16.,64.])[:,None]/(k*t))
    gx = np.sort(np.clip(np.concatenate(spans),-1.,1.),axis=0).T
    dx = np.diff(gx,axis=1)[:,:,None]/2
    live = dx>0.
    X = (gx[:,1:]+gx[:,:-1])[:,:,None]/2 + dx*fm
    X = np.where(live,X,0.)
    Jx=0.
    if mode=='rectangle':
        s0=seg(k*t*X,k*t*X/2)
        qu = k*mu*r
        J0=t*s0*(-np.expm1(-1j*qu))
        Jz=r*seg(-qu,-qu/2)*np.expm1(1j*k*t*X)
    else:
        a=k*t*(X-v*mu)
        b=a-k*(mu*r+dr*rx)
        A=seg(a,a/2)
        B=seg(b,b/2)
        J0=t*(A-B)
        joined=seg(b-a,(a+b)/2)-B
        Jz=v*J0+r*joined
        Jx=rx*joined
    err=np.max(abs(k*(X*J0-mu*Jz-dr*Jx)))/(1+k*t)
    assert err < 5e-14
    ls=1-X*np.arctanh(X)
    li=np.pi*X/2
    ts=.5*(X*X+X*(1-X*X)*np.arctanh(X))
    ti=np.pi*X*(1-X*X)/4
    L=abs(J0)**2/((k*k+ls)**2+li*li)
    jt=abs(np.sqrt(1-mu*mu)*Jz-mu*np.cos(az)*Jx)**2+abs(np.sin(az)*Jx)**2
    T=jt*.5*(1-X*X)/((k*k*(1-X*X)+ts)**2+ti*ti)
    norm=dx*fw*aw[:,None,None]/2
    return np.array([np.sum(norm*L),np.sum(norm*T)])*k*k/(2*np.pi*t)


def fold(mode,t,r,v,nr,na,nf,kmax,low=3e-4,rx=0.,az=0.):
    z,w=gauss(nr)
    k=np.r_[0.,np.geomspace(low,kmax,34)]
    q=(k[1:]+k[:-1])[:,None]/2+np.diff(k)[:,None]*z/2
    qw=(np.diff(k)[:,None]*w/2).flatten()
    rows=np.array([radial_loop(p,mode,t,r,v,na,nf,rx=rx,az=az) for p in q.flatten()])
    return qw@rows


def run(level):
    nr,na,nf=(3,24,12) if level==0 else (5,38,20) if level==1 else (6,58,32)
    final_r=.601034457155713
    # Standalone radius scan in units m_D; t is in 1/m_D.
    rows=[]
    for kt in (8.,20.):
        for mode,t,r,v in [('rectangle',32.,final_r,0.),('rectangle',160.,final_r,0.),
                           ('triangle',32.,final_r,.97),('triangle',160.,final_r,.994),
                           ('triangle',500.,final_r,.998),
                           ('triangle',996.,final_r,.9982),
                           ('triangle',32.,0.,.97)]:
            ans=fold(mode,t,r,v,nr,na,nf,kt)
            assert min(ans)>=0.
            if r==0.: assert max(ans)==0.
            ref=static_shape(r,False,120)[0] if mode=='rectangle' else None
            local=weights(v,r)[1:].tolist() if mode=='triangle' else None
            rows.append(dict(mode=mode,t=t,final_r=r,v=v,kmax=kt,
                             finite_L=ans[0],finite_T=ans[1],
                             longitudinal_rectangular_static=ref,
                             doppler_local_triangle=local))
    return dict(scope='Classical symmetrized spacelike HTL only; closed straight-path contours; '
                      'no timelike pole, hard transfer completion, or kinetic rate.',level=level,rows=rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--level',type=int,default=0)
    print(json.dumps(run(p.parse_args().level),indent=2,allow_nan=False))
