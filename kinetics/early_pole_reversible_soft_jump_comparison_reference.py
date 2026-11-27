"""Stationary projected radial jump vs matched diffusion, with a soft cutoff."""
import argparse
import json
import numpy as np
from scipy.special import roots_legendre
from scipy.sparse import diags,eye,kron
from scipy.sparse.linalg import expm_multiply,spsolve

from early_pole_kick_characteristic_reference import sphere
from early_pole_matched_wilson_geometry_reference import inputs
import early_pole_finite_resonance_aperture_reference as q


def rates(g,cut,h):
    d=inputs(g)
    a,mD=d['angular_halfspan'],g*np.sqrt(2/3)
    W,kz,kp=sphere(cut,10,d['vP'])
    # Symmetric transfer. Our positive-axis kernel accounts for both signs.
    W=W*g*g/(8*np.pi*a)
    dl=np.broadcast_to(mD*kz/h,W.shape)
    low=np.floor(dl).astype(int)
    part=(dl*dl-low*low)/(2*low+1)
    bins=np.bincount(low.ravel(),weights=(W*(1-part)).ravel(),minlength=int(low.max())+2)
    bins+=np.bincount((low+1).ravel(),weights=(W*part).ravel(),minlength=bins.size)
    moment=np.sum(bins*(np.arange(bins.size)*h)**2)
    exact=np.sum(W*(kz*mD)**2)
    assert abs(moment/exact-1)<2e-13
    trans=np.sum(W*(kp*mD)**2)*2
    return bins,exact/h**2,trans


def calc(g,cut,n,domain,source_width,nang):
    U0,EL,EH,_,_=q.intervals(np.hypot(q.X,g/2),q.MS,np.hypot(q.RB*q.X,g/2),q.Q)
    M,B=np.hypot(q.X,g/2),np.hypot(q.RB*q.X,g/2)
    center=np.sqrt(((EL+EH)/2-M)*((EL+EH)/2+M))
    d=inputs(g)
    ang,vw=roots_legendre(nang)
    vw=vw/2
    h=domain/n
    p=center-domain/2+h*(np.arange(n)+.5)
    ep=np.hypot(p,M)
    R=p*p*np.exp(-ep); R/=R.sum()
    z=np.sqrt(np.outer(R,vw)).ravel()
    u=(np.sqrt(np.outer(R,vw))*np.exp(-.25*((p-center)/source_width)**2)[:,None]).ravel()
    u/=np.linalg.norm(u)
    k=p[:,None]
    eb=np.sqrt(k*k+q.Q*q.Q-2*k*q.Q*ang+B*B)
    defect=((M*M-B*B-q.Q*q.Q+2*k*q.Q*ang)/(ep[:,None]+eb)-U0).ravel()/d['angular_halfspan']
    jump,D,trans=rates(g,cut,h)
    # Match the angular dipole rate; isotropic redraw is illustrative, not the HTL angular operator.
    nu=trans/(2*center*center)
    cut_rho=[]
    for name,taps in [('ballistic',np.array([0.])),('nonlocal',jump),('diffusive',np.array([0.,D]))]:
        offsets=[0]
        links=[np.zeros(n)]
        for j,kick in enumerate(taps[1:],1):
            if j>=n:break
            links[0][:-j]-=kick*np.sqrt(R[j:]/R[:-j])
            links[0][j:]-=kick*np.sqrt(R[:-j]/R[j:])
            offsets.extend((-j,j))
            links.extend((np.full(n-j,kick),np.full(n-j,kick)))
        F=kron(diags(links,offsets,format='csc'),eye(nang),format='csc')
        assert np.max(abs(F@z))<3e-12
        X=F-1j*diags(defect,format='csc')
        fx=expm_multiply(X,u,start=0.,stop=4.,num=9,traceA=X.diagonal().sum())
        C=fx@u
        assert abs(C[0]-1.)<1e-12 and max(abs(C))<1+1e-12
        D0=.35*eye(len(z),format='csc')-F
        # Stationary line with a common small readout regulator; check positive quadratic identity.
        spec=[]
        for w in (-2.,0.,2.):
            y=spsolve(D0+1j*diags(defect-w,format='csc'),u)
            v=np.vdot(u,y).real
            assert v>=0. and abs(v-np.vdot(y,D0@y).real)<2e-9
            spec.append(v)
        cut_rho.append(dict(case=name,abs_correlation=abs(C).tolist(),
                            correlation_real=C.real.tolist(),correlation_imag=C.imag.tolist(),
                            line_at_w_minus2_zero_plus2=spec))
    return dict(g=g,cut=cut,n=n,domain=domain,source_width=source_width,nang=nang,
                radial_bulk_diffusion_jump=D,angular_redraw_rate=nu,jump_bins_ignoring_zero=len(jump)-1,rows=cut_rho)


def run(level):
    n,nang=(49,22) if level==0 else (83,32) if level==1 else (119,40)
    return dict(scope='Prescribed stationary soft kick generator. Equal bulk second moment, '
                      'no exact recoil, vertex or thermal production rate.',
                time_in_inverse_angular_halfspans=np.linspace(0.,4.,9).tolist(),
                rows=[calc(g,c,n,.28,w,nang) for g,c,w in
                      [(.005,10.,.008),(.005,20.,.008),(.005,20.,.014),(.003,20.,.008)]])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--level',type=int,default=0)
    print(json.dumps(run(p.parse_args().level),indent=2,allow_nan=False))
