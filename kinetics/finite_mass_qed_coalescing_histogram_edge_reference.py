"""Cumulative double-root scaling at actual archive edges, diagnostic widths only."""
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import roots_legendre

from early_pole_hard_yukawa_photoproduction_reference import I,G,dot,slash,spinors
from finite_mass_qed_frozen_transverse_plasmon_gate_reference import Q
from finite_mass_qed_nonnull_photon_ward_poles_reference import setup


def cumulative(y,h):
    if h==0:
        logpart=2*y/(1+y*y)
    else:
        logpart=np.log1p(h*(2*y-h)/(1+(y-h)**2))/h
    return np.pi+np.arctan(y)+np.arctan(y-h)+logpart


def weight(y,w,b):
    return (w[0]-w[1])*cumulative(y,b)+w[1]*2*np.pi


def actual(d,g,h,eta,e):
    k,c=.0008,.25
    Es=np.hypot(Q,d['checkpoints'][-1]['scalar_pole_over_T'])
    omega=(h*e/2+Q*k*c)/Es
    M,m,q,K,T,geo=setup(d,k,omega,c)
    nm,idx=('fermion',32) if g=='fine' else ('parent',48)
    ref=np.hypot(m if g=='fine' else M,d[f'endpoint_{nm}_momentum_edges'][idx])
    if g=='fine':ref+=T[0]
    V=np.linalg.norm(T[1:]);ax=T[1:]/V;cross=np.array([ax[2],0.,-ax[0]])
    C=(M*M+dot(T,T)-m*m)/2
    signed_K=dot(K,K)
    def dd(E,p):return (signed_K+2*dot(geo(E,p)[0],K))/e
    ph=brentq(lambda p:dd(ref,p)-eta,0,np.pi,xtol=2e-15)
    Z=slash(K)
    def at(x):
        ep=brentq(lambda z:dd(z,ph)-x,ref-.003,ref+.003,xtol=3e-16)
        p,b=geo(ep,ph)
        v=(ep*T[0]-C)/V
        u=np.sqrt((ep-M)*(ep+M)-v*v)
        deriv=ax[2]*T[0]/V+cross[2]*np.cos(ph)*(ep-v*T[0]/V)/u
        jac=1/(2*(omega-k*deriv))
        if not jac>0:raise ArithmeticError('Degenerate first-root chart')
        spin=spinors(p,M);bar=spinors(b,m).conj().T@G[0]
        A=slash(p)+Z+M*I;B=slash(b)-Z+m*I
        left=np.ravel([bar@A@slash(z)@spin for z in np.eye(4)[1:3]])
        right=np.ravel([bar@slash(z)@B@spin for z in np.eye(4)[1:3]])
        i,j=[np.searchsorted(d[f'endpoint_{name}_momentum_edges'],np.linalg.norm(X[1:]))-1
             for name,X in (('parent',p),('fermion',b))]
        fp=d['endpoint_parent_occupations'][i];fb=d['endpoint_fermion_occupations'][j]
        return jac,left,right,fp*(1-fb)
    draw=2*dot(q,K)/e
    if abs(draw-h)>3e-11 or abs(dd(ref,ph)-eta)>3e-6:
        raise ArithmeticError('Requested joint root/edge scaling not recovered')
    J,A,B,_=at(0.)
    norm=J*np.sum(abs(A)**2)
    diff=np.sum(abs(A-B)**2)/np.sum(abs(A)**2)
    lower,upper=at(eta-1.)[-1],at(eta+1.)[-1]
    if lower<=0 or upper<=0 or lower==upper:raise ArithmeticError('Real histogram edge missing')
    mean=(lower+upper)/2
    W=np.array([lower,upper])/mean
    analytic=weight(eta,W,h)
    convert=h*h/(h*h+4)
    xs,ws=roots_legendre(50)
    guards=np.unique(np.sort(np.r_[-320,320,eta, np.linspace(-8, h+8,12),-64,-24,h+24,h+64]))
    vals=np.zeros(2)
    for a,b in zip(guards[:-1],guards[1:],strict=True):
        for x,wj in zip((a+b)/2+(b-a)*xs/2,(b-a)*ws/2,strict=True):
            J,A,B,f=at(x)
            if f!=[lower,upper][int(x>eta)]:
                raise ArithmeticError('Numerical integral crossed a second histogram edge')
            m0=A/(x+1j)
            vals+=wj*(J*f/(norm*mean))*np.array([np.sum(abs(m0+B/(draw-x+1j))**2),
                                                  np.sum(abs(m0+B/(draw-x-1j))**2)])
    pred=np.array([analytic,analytic*convert])
    E=float(max(abs(vals-pred)))
    if E>7e-4:
        raise ArithmeticError(('Joint scaling at actual edge',g,h,eta,e,vals,pred,diff))
    return dict(grid=g,separation_in_widths=h,edge_displacement_in_widths=eta,width=e,
                physical_k=k,physical_omega=omega,edge_parent_energy=ref,azimuth=ph,
                normalized_step_weights=W.tolist(),relative_residue_difference_squared=diff,
                actual_same_and_opposite_imag_sign=vals.tolist(),joint_limit=pred.tolist(),
                maximum_joint_deviation_absolute=E)


def run():
    out=[];d0=0.
    for h in (0.,.5,2.,6.,16.):
        for eta in (-2.,0.,3.):
            v=quad(lambda x:(h*h+4)/((x*x+1)*((h-x)**2+1)),-np.inf,eta,epsabs=1e-12,epsrel=1e-12)[0]
            d0=max(d0,abs(v-cumulative(eta,h)))
    if d0>5e-12:raise ArithmeticError('Normalizing analytic two-root primitive')
    path=Path(__file__).parent/'results'
    for g in ('fine','third'):
        with gzip.open(path/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as f:d=json.load(f)
        for h in (0.,.5,2.,6.,16.):
            for x in (-2.,0.,3.):
                for e in (2e-10,4e-11):
                    out.append(actual(d,g,h,x,e))
    return dict(scope='Joint coalescence/step diagnostics with declared imaginary signs only.',
                exact_cumulative_leading_primitive_defect=d0,data=out)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
