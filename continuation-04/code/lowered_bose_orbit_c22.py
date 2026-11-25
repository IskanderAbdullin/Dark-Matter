"""Initial orbit-averaged C22 relative rate and radial periods for occupied s>0."""
import argparse,json,numpy as np
from scipy.optimize import brentq
from lowered_antonov_hessian import Halo
from lowered_bose_signed_source import nodes
from lowered_bose_pointwise_c22 import c22
from lowered_bose_exterior_source import R,V,GYR,TIME

def orbit(H,s,z,n=24,nr=24):
    rs=brentq(lambda r:float(H.psi(r))-s,0,H.rt,xtol=1e-13)
    rc=brentq(lambda r:2*(float(H.psi(r))-s)-float(H.sol.sol(max(r,H.eps))[1])/max(r,H.eps),H.eps,rs,xtol=1e-13)
    lm=2*rc**2*(float(H.psi(rc))-s);ell=lm*z
    f=lambda r:2*r*r*(float(H.psi(r))-s)-ell
    rp=brentq(f,0,rc,xtol=1e-13) if z>0 else 0.;ra=brentq(f,rc,rs,xtol=1e-13) if z>0 else rs
    x,w=nodes(nr);th=np.pi*x/2;r=rp+(ra-rp)*np.sin(th)**2;dr=(ra-rp)*np.pi*np.sin(th)*np.cos(th)*w
    ps=H.psi(r);q=np.sqrt(np.maximum(ps-s,0));rad=np.sqrt(2*(ps-s)-ell/r**2);wt=dr/rad
    raw=np.array([c22(H,p,qi,n) for p,qi in zip(ps,q)])
    rho=V**2/(4*np.pi*6.6743e-11*R**2);pref=rho*.125*np.sqrt(2)*V*GYR/H.Ic
    rate=pref*raw/float(H.F(s));mean=wt@rate/wt.sum();period=2*TIME*wt.sum()
    return dict(s=s,gap=H.eta-s,L2frac=z,rp=rp,ra=ra,period_Gyr=period,rate_Gyr=mean,local_min=float(rate.min()),local_max=float(rate.max()),change_per_period=period*mean)

def run(H,n,nr):
    rows=[orbit(H,s,z,n,nr) for s in [.0001,.001,.01,.1,.5,1.,2.,3.,3.5,3.9,3.99,3.999,3.9999] if s<H.psi0 for z in [.01,.5,.99]]
    return dict(model=H.model,n=n,nr=nr,rows=rows)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--model',default='linear');p.add_argument('--n',type=int,default=24);p.add_argument('--nr',type=int,default=24);a=p.parse_args();print(json.dumps(run(Halo(model=a.model),a.n,a.nr)))
