"""Conditional isolated homogeneous cell equilibrium criterion; not halo condensation."""
import argparse,json,numpy as np
from scipy.special import gamma,zeta
from scipy.optimize import brentq
from lowered_antonov_hessian import Halo
from lowered_bose_signed_source import nodes

def run(H,n=256):
    x,w=nodes(n);cn=gamma(1.5)*zeta(1.5)/2;ce=gamma(2.5)*zeta(2.5)/2
    def values(p):
        p=np.asarray(p);q=np.sqrt(p[...,None])*x;f=H.F(p[...,None]-q*q);N=H.I(p);E=p**2.5*np.sum(w*x**4*f,axis=-1);beta=(ce/E)**.4;Nc=cn/beta**1.5
        return N,E,Nc,beta
    def ratio(p):
        N,E,Nc,beta=values(p);return N/Nc
    if ratio(H.psi0)<=1:return dict(model=H.model,delta=H.delta,n=n,center_ratio=float(ratio(H.psi0)),supercritical=False)
    p=brentq(lambda p:ratio(p)-1,1e-4,H.psi0);r=brentq(lambda r:float(H.psi(r))-p,0,H.rt);rr=r*x;N,E,Nc,beta=values(H.psi(rr));ex=r*np.sum(w*rr*rr*np.maximum(N-Nc,0))/(H.Ic*H.Qt)
    return dict(model=H.model,delta=H.delta,n=n,center_ratio=float(ratio(H.psi0)),supercritical=True,threshold_psi=p,radius_R=r,mass_fraction=float(H.sol.sol(r)[1]/H.Qt),conditional_excess_fraction=float(ex),center_excess_fraction=float(1-1/ratio(H.psi0)))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--model',default='linear');p.add_argument('--delta',type=float,default=0);p.add_argument('--n',type=int,default=256);a=p.parse_args();print(json.dumps(run(Halo(model=a.model,delta=a.delta),a.n)))
