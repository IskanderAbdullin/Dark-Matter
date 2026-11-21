"""Numerical excess-mass mixing diagnostic for a finite critical lowered-Bose halo.
Requires numpy scipy and sibling critical_lowered.py.
Uses 6-D phase-space superlevel volumes, unlike comparisons of only density or momentum.
No collision solver or actual cosmological initial conditions; initial patch is hypothetical.
"""
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq
from scipy.special import zeta
from critical_lowered import solve

C0=4*np.sqrt(3)*np.pi**3
Cm=3*np.sqrt(3)*zeta(3)*np.sqrt(np.pi/2)

def legint(fun,upper,n=48):
    xx,ww=leggauss(n);u=(xx+1)*upper/2
    return upper/2*np.dot(ww,fun(u))

def vol_ratio(L,eta,rt,sol,B):
    if L<3e-6:return 1+B*np.sqrt(L)
    def y(r):return eta-sol.sol(r)[0]
    rF=brentq(lambda r:float(y(r)-L),1e-5,rt,xtol=1e-12)
    xx,ww=leggauss(72);rs=(xx+1)*rF/2
    integral=rF/2*np.dot(ww,rs*rs*np.maximum(L-y(rs),0)**1.5)
    return (16*np.pi**2/3)*2**1.5*integral/(C0*L**3)

def excess(F,eta,rt,sol,Ic,Q,n=48):
    nu=4*Ic/np.sqrt(np.pi)
    s=nu*Q/Cm
    B=32/(35*Ic)
    c=1/np.expm1(eta)
    Lin=np.log1p(1/F);Lfin=np.log1p(1/(F+c))
    kernel=lambda x:x**3*np.exp(x)/np.expm1(x)**2
    Di=s*legint(kernel,Lin,n)
    xx,ww=leggauss(n);ls=(xx+1)*Lfin/2
    v=np.array([vol_ratio(float(x),eta,rt,sol,B) for x in ls])
    Df=Lfin/2*np.dot(ww,kernel(ls)*v)
    return dict(F=F,Dfinal_Dinitial=float(Df/Di),
                initial_leading_phase_volume_over_final=float(s))

def run():
    out=[]
    for eta in (.5,1,1.35,1.4,1.5,2,3,4,6):
        Ic,rt,Q,sol,_=solve(eta,nq=192)
        nu=4*Ic/np.sqrt(np.pi)
        thresholds=[excess(F,eta,rt,sol,Ic,Q)
                    for F in (.0001,.001,.01,.1,1,10,100,1000,10000)]
        out.append(dict(eta=eta,nu=nu,Q=Q,
                        highF_mass_margin=nu*Q/Cm,thresholds=thresholds))
    e4=next(e for e in out if e['eta']==4)
    val=next(t['Dfinal_Dinitial'] for t in e4['thresholds'] if t['F']==.001)
    assert 1.015<val<1.017
    assert e4['highF_mass_margin']>8
    assert next(e for e in out if e['eta']==1)['highF_mass_margin']<1
    return dict(model='closed thermal real-scalar patch -> exactly critical lowered Bose, no escape',
       highF_necessary='nu(eta)*Q(eta) >= 3sqrt(3)zeta(3)sqrt(pi/2)',
       general_lowF='Dfinal> Dinitial for sufficiently small F: final compact phase support vs initial unbounded thermal momentum tail, both with same particle number',
       tested=3,rows=out,
       caution='All limits are conditional on initial thermal spectrum, mass closure, scalar occupation, and absence of collisions/coherence. No real galaxy excluded.')
if __name__=='__main__': print(json.dumps(run(),indent=2))
