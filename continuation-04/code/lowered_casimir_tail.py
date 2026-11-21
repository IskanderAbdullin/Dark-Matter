"""Precision high-F six-dimensional phase-volume of an exactly critical lowered BE halo.
Independent small-r high-accuracy potential integration; not a kinetic or galaxy simulation.
Requires scipy numpy and sibling critical_lowered.py.
"""
import json
import numpy as np
from scipy.integrate import solve_ivp,quad
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from critical_lowered import I

def run(eta=4.,nq=300):
    Ic=I(eta,eta,nq);C=np.pi/(2*Ic)
    ys=np.r_[0,np.geomspace(1e-13,.2,1900)]
    xx=I(eta-ys,eta,nq)/Ic
    xx[:5]=1-C*np.sqrt(ys[:5])
    eos=PchipInterpolator(ys,xx)
    def density(y):
        return max(0,1-C*np.sqrt(y)) if y<1e-10 else float(eos(y))
    A=np.pi/(2*np.sqrt(6)*Ic)
    eps=1e-6
    def ode(r,w):return[w[1]/r**2,r*r*density(max(w[0],0))]
    initial=[eps**2/6-A*eps**3/12,eps**3/3-A*eps**4/4]
    res=solve_ivp(ode,[eps,2],initial,method="DOP853",rtol=1e-12,
                  atol=1e-15,max_step=.0005,dense_output=True)
    assert res.success
    B=32/(35*Ic)
    rows=[]
    for F in (100,1000,10000,100000,1000000):
        c=1/np.expm1(eta)
        L=np.log1p(1/(F+c))
        R=brentq(lambda r:float(res.sol(r)[0]-L),eps,2,xtol=1e-12)
        volume=(16*np.pi**2/3)*2**1.5*quad(
          lambda r:r*r*max(L-res.sol(r)[0],0)**1.5,0,R,
          epsabs=1e-17,epsrel=1e-10)[0]
        leading=4*np.sqrt(3)*np.pi**3*L**3
        ratio=volume/leading
        rows.append(dict(F=F,exact_ratio=ratio,
            predicted_first_order=1+B*np.sqrt(L),
            error=ratio-1-B*np.sqrt(L),rF=R))
    assert abs(rows[-2]['error'])<1e-5
    assert abs(rows[-1]['error'])<2e-6
    assert abs(rows[-1]['exact_ratio']-1)<.002
    return dict(eta=eta,Ic=Ic,B=B,checks_passed=3,rows=rows,
                warning="F->infinity static 6D phase-volume asymptotic, not an observed velocity distribution.")
if __name__=="__main__": print(json.dumps(run(),indent=2))
