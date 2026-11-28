"""Draft stationary line from an equilibrium reversible momentum-jump chain."""
import json

import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre
from scipy.sparse import csc_matrix, diags, eye
from scipy.sparse.linalg import spsolve

import early_pole_finite_resonance_aperture_reference as q


def model(radial, angular, size, nang=8, extent=.15, source_width=None):
    ga, aw=roots_legendre(nang)
    aw=aw/2
    p=q.PC-extent/2+(.5+np.arange(size))*extent/size
    ep=np.hypot(p,q.MP)
    R=p*p*np.exp(-ep/.8)
    R/=R.sum()
    u=np.sqrt(np.outer(R,aw)).ravel()
    # Preserve the imposed diffusion strength as the spacing changes.
    jump=radial*(size/9*.15/extent)**2
    L=np.zeros((size,size))
    for i in range(size-1):
        f=jump*np.sqrt(R[i]*R[i+1])
        L[i,i]-=f/R[i]
        L[i+1,i+1]-=f/R[i+1]
        L[i+1,i]=L[i,i+1]=jump
    za=np.sqrt(aw)
    A=angular*(np.outer(za,za)-np.eye(aw.size))
    F=np.kron(L,np.eye(aw.size))+np.kron(np.eye(p.size),A)
    assert max(abs(F@u))<2e-10
    k=p[:,None]
    B=np.sqrt(k*k+q.Q*q.Q-2*k*q.Q*ga+q.MB*q.MB)
    d=(q.MP*q.MP-q.MB*q.MB-q.Q*q.Q+2*k*q.Q*ga)/(ep[:,None]+B)-q.ES
    d=d.ravel()/(np.diff(q.defects(q.PC))[0]/2)
    u=(np.sqrt(np.outer(R,aw))*(np.exp(-.25*((p-q.PC)/source_width)**2)[:,None] if source_width else 1.)).ravel()
    u/=np.linalg.norm(u)
    D=.002*eye(len(d),format='csc')-csc_matrix(F)
    Dx=D+1j*diags(d,format='csc')
    Id=eye(len(d),format='csc')
    def line(x):
        y=spsolve(Dx-1j*x*Id,u)
        v=np.vdot(u,y).real
        direct=np.vdot(y,D@y).real
        assert abs(v-direct)<5e-11
        return v
    a=quad(line,-100,100,epsabs=4e-7,points=np.linspace(-15,15,43),limit=400)[0]
    a+=quad(line,-np.inf,-100,epsabs=1e-8)[0]+quad(line,100,np.inf,epsabs=1e-8)[0]
    assert abs(a/np.pi-1)<1e-5
    return dict(radial=radial,actual_per_link_jump=jump,angular=angular,size=size,angular_nodes=nang,domain_extent=extent,source_width=source_width,zero=line(0.),
                width_proxy=(a/np.pi)/line(0.),normalized_area=a/np.pi)


if __name__=='__main__':
    print(json.dumps(dict(scope='Stationary reversible momentum redraw toy; no kinetic calibration.',
            rows=[model(r,2.,n) for n in (9,15,27,39) for r in (.1,1.,10.,100.)]+
                 [model(r,2.,39,nang=16) for r in (.1,1.,10.,100.)]+
                 [model(r,2.,n,nang=12,extent=x,source_width=.008) for n,x in ((31,.15),(51,.23),(61,.30))
                  for r in (.1,1.,10.,100.)]),indent=2,allow_nan=False))
