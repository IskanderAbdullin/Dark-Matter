"""Independent fixed-(s,L²) radial covariance susceptibility, weak FE basis."""
import argparse,json,numpy as np
from scipy.optimize import brentq
from scipy.linalg import eigh
from lowered_antonov_hessian import Halo
from lowered_bose_signed_source import nodes

def run(H,nf=32,ns=64,nl=24,nr=64):
    grid=H.rt*np.linspace(0,1,nf+1)**2
    A=np.zeros((nf,nf));B=np.zeros_like(A);shell=np.zeros_like(A)
    for i,(a,b) in enumerate(zip(grid[:-1],grid[1:])):
        val=(b**3-a**3)/(3*(b-a)**2)
        A[i,i]+=val
        if i+1<nf:A[i+1,i+1]+=val;A[i,i+1]-=val;A[i+1,i]-=val
    def basis(r):
        k=np.clip(np.searchsorted(grid,r,side='right')-1,0,nf-1);t=(r-grid[k])/(grid[k+1]-grid[k]);X=np.zeros((len(r),nf+1));j=np.arange(len(r));X[j,k]=1-t;X[j,k+1]=t;return X[:,:nf]
    def cov(r,w):
        X=basis(r);m=w@X;return (X.T*w)@X-np.outer(m,m)/w.sum()
    u,wu=nodes(ns);z,wz=nodes(nl);v,wv=nodes(nr);max_phase_error=0;sum_phase_error=0;sum_phase=0;min_rad=1e99
    for ui,wi in zip(u,wu):
        s=H.psi0*(1-ui**2);ds=2*H.psi0*ui*wi
        b=1/np.expm1(H.eta-s);fp=b*(1+b)
        if H.model=='quadratic':fp-=H.d
        rs=brentq(lambda r:float(H.psi(r))-s,0,H.rt,xtol=1e-13)
        rc=brentq(lambda r:2*(float(H.psi(r))-s)-float(H.sol.sol(max(r,H.eps))[1])/max(r,H.eps),H.eps,rs,xtol=1e-13) if rs>2*H.eps else rs/np.sqrt(2)
        lm=2*rc*rc*(float(H.psi(rc))-s)
        cuts0=np.r_[0,grid[(grid>0)&(grid<rs)],rs];r=(cuts0[:-1,None]+np.diff(cuts0)[:,None]*v).ravel();wr=(np.diff(cuts0)[:,None]*wv).ravel();w=wr*r*r*np.sqrt(np.maximum(H.psi(r)-s,0))/2
        shell+=ds*fp/H.Ic*cov(r,w);p0=w.sum();p1=0
        for zi,wli in zip(z,wz):
            ell=lm*zi
            f=lambda r:2*r*r*(float(H.psi(r))-s)-ell
            rp=brentq(f,0,rc,xtol=1e-13);ra=brentq(f,rc,rs,xtol=1e-13)
            inside=grid[(grid>rp)&(grid<ra)];cuts=np.r_[0,np.arcsin(np.sqrt((inside-rp)/(ra-rp))),np.pi/2]
            th=(cuts[:-1,None]+np.diff(cuts)[:,None]*v).ravel();wt=(np.diff(cuts)[:,None]*wv).ravel()
            rr=rp+(ra-rp)*np.sin(th)**2;dr=2*(ra-rp)*np.sin(th)*np.cos(th)*wt
            rad=H.psi(rr)-s-ell/(2*rr*rr);min_rad=min(min_rad,rad.min())
            if np.any(rad<=0):raise ValueError((s,zi,rad.min(),rp,ra))
            w=dr/(8*np.sqrt(rad))*lm*wli
            B+=ds*fp/H.Ic*cov(rr,w);p1+=w.sum()
        err=abs(p1/p0-1);max_phase_error=max(max_phase_error,err);sum_phase_error+=ds*fp*abs(p1-p0);sum_phase+=ds*fp*p0
    eb=eigh(B,A,eigvals_only=True);es=eigh(shell,A,eigvals_only=True);ed=eigh(shell-B,A,eigvals_only=True)
    return dict(model=H.model,delta=H.delta,nf=nf,ns=ns,nl=nl,nr=nr,orbit_max=float(eb[-1]),shell_max=float(es[-1]),orbit_min=float(eb[0]),shell_minus_orbit_min=float(ed[0]),phase_max_rel=float(max_phase_error),phase_weighted_rel=float(sum_phase_error/sum_phase),min_radicand=float(min_rad),orbit_top5=eb[-5:].tolist(),shell_top5=es[-5:].tolist())
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--model',default='linear');p.add_argument('--delta',type=float,default=0);p.add_argument('--nf',type=int,default=32);p.add_argument('--ns',type=int,default=64);p.add_argument('--nl',type=int,default=24);p.add_argument('--nr',type=int,default=16);a=p.parse_args();print(json.dumps(run(Halo(model=a.model,delta=a.delta),a.nf,a.ns,a.nl,a.nr)))
