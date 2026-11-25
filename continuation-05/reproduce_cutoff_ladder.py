"""Positive Bose Taylor cutoffs and radial Poisson profiles.
Python 3 + numpy, scipy; no network, GitHub Actions or formation/stability claims.
Run: python reproduce_cutoff_ladder.py
"""
import json
import numpy as np
from math import factorial,pi,log,exp
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from numpy.polynomial.legendre import leggauss
N,W=leggauss(160);U=(N+1)/2;W=W/2

def occupation(t,eta,k):
    t=np.asarray(t)
    n=np.arange(1,100,dtype=float)
    c=np.array([(n**j*np.exp(-eta*n)).sum() for j in range(k+11)])
    out=np.zeros_like(t,dtype=float)
    near=t<.16
    if np.any(near):
        a=t[near]
        out[near]=sum(c[j]*a**j/factorial(j) for j in range(k,k+11))
    if np.any(~near):
        a=t[~near]
        out[~near]=1/np.expm1(eta-a)-sum(c[j]*a**j/factorial(j) for j in range(k))
    return out

def I(p,eta,k):
    p=np.asarray(p,dtype=float)
    return p**1.5*np.sum(W*U**2*occupation(p[...,None]*(1-U*U),eta,k),axis=-1)

def halo(k,psi0=4.,eta=4.,max_radius=1e10):
    x=np.unique(np.r_[np.geomspace(1e-12,psi0,550),
             np.linspace(.0001,psi0,450),
             psi0-np.geomspace(1e-10,psi0,280)])
    x=x[x>0];y=I(x,eta,k);Ic=float(y[-1])
    fun=PchipInterpolator(np.log(x),np.log(y/Ic))
    xmin=x[0];a=1e-4
    def rho(p):
        if p<=0:return 0.
        if p<xmin:return float(np.exp(fun(log(xmin)))*(p/xmin)**(k+1.5))
        return float(np.exp(fun(log(min(p,psi0)))))
    def rhs(t,z):
        psi,v=z
        return [-v,exp(2*t)*rho(psi)-v]
    def edge(t,z):return z[0]
    edge.terminal=True;edge.direction=-1
    sol=solve_ivp(rhs,[log(a),log(max_radius)],
        [psi0-a*a/6,a*a/3],events=edge,method='DOP853',
        first_step=.01,max_step=.12,rtol=4e-9,atol=2e-11,dense_output=True)
    hit=len(sol.t_events[0])>0
    te=float(sol.t_events[0][0]) if hit else log(max_radius)
    psi,v=sol.sol(te);r=exp(te)
    out={'k':k,'psi0':psi0,'finite_edge_found':hit,
         'radius_over_R_or_limit':r,'Q_at_end':float(v*r),
         'nu':4*Ic/np.sqrt(pi)}
    if not hit:
        rr=np.logspace(7,10,25)
        z=sol.sol(np.log(rr));m=rr*z[1]
        out['mass_log_slope']=float(np.polyfit(np.log(rr),np.log(m),1)[0])
    return out

def main():
    tests=[]
    for k in range(1,7):
        xx=np.geomspace(1e-7,3.99,160)
        ff=occupation(xx,4,k)
        assert np.all(ff>0) and np.all(np.diff(ff)>0)
        px=np.array([.005,.01,.02,.04]);ii=I(px,4,k)
        slope=np.polyfit(np.log(px),np.log(ii),1)[0]
        assert abs(slope-(k+1.5))<.04
        tests.append({'k':k,'edge_density_exponent':float(slope)})
    out=[halo(k,psi0=4.) for k in range(1,7)]
    out += [halo(2,psi0=3.8),halo(3,psi0=3.8)]
    assert abs(out[0]['radius_over_R_or_limit']-40.2988)<.02
    assert abs(out[1]['radius_over_R_or_limit']-130.178)<.1
    assert out[2]['finite_edge_found'] and out[2]['radius_over_R_or_limit']>1e8
    for k in (4,5,6):
        d=out[k-1];assert not d['finite_edge_found']
        n=k+1.5;predict=(n-3)/(n-1)
        assert abs(d['mass_log_slope']-predict)<.08
    return {'positive_monotone_edge_tests':tests,'profiles':out,
            'analytic_edge_exponent':'k+3/2',
            'limitations':'No global proof of no radius for k>=4 and no nonlinear stability or physical formation theory'}
if __name__=='__main__':print(json.dumps(main(),indent=2))
