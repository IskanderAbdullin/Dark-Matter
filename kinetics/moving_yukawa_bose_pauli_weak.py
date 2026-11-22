"""Moving on-shell Dirac χ1<->χ2+real φ scalar Yukawa Pauli/Bose weak source.

Fully boosted two-body daughters, exact physical threshold and scalar trace;
conditional one-way+inverse quantum source for prescribed f1,f2,bphi.
Numerical synthetic phase-space control only, NOT strong nonlinear cosmological
χ/φ dynamics, physical Γelastic/FRW, BEC or 86eV DM conclusion.
"""
from __future__ import annotations
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import expit


def kallen(m1,m2,mphi):
    if not all(math.isfinite(x) and x>=0 for x in (m1,m2,mphi)) or m1<=0:
        raise ValueError('finite nonnegative masses, positive mother')
    if m1<=m2+mphi:
        return 0.
    return (m1*m1-(m2+mphi)**2)*(m1*m1-(m2-mphi)**2)


def vacuum_width(m1,m2,mphi,y):
    if not math.isfinite(y) or y<0:raise ValueError('nonnegative coupling')
    lam=kallen(m1,m2,mphi)
    return y*y*((m1+m2)**2-mphi*mphi)*math.sqrt(lam)/(16*math.pi*m1**3) if lam>0 else 0.


def daughter_energies(mother_p,c,m1,m2,mphi):
    if not(math.isfinite(mother_p) and mother_p>=0 and
           np.all(np.isfinite(c)) and np.all(np.asarray(c)>=-1) and
           np.all(np.asarray(c)<=1)):
        raise ValueError('physical mother and angular coordinates')
    lam=kallen(m1,m2,mphi)
    if lam<=0: raise ValueError('closed on-shell decay')
    P=math.hypot(m1,mother_p)
    e2star=(m1*m1+m2*m2-mphi*mphi)/(2*m1)
    phistar=(m1*m1-m2*m2+mphi*mphi)/(2*m1)
    pstar=math.sqrt(lam)/(2*m1)
    boost=mother_p/m1
    gamma=P/m1
    e2=gamma*e2star+boost*pstar*np.asarray(c)
    eph=gamma*phistar-boost*pstar*np.asarray(c)
    return e2,eph


def event(p,c,m1,m2,mphi,f1,f2,bphi):
    e1=math.hypot(p,m1)
    e2,eph=daughter_energies(p,c,m1,m2,mphi)
    a=float(f1(e1));d=np.asarray(f2(e2));b=np.asarray(bphi(eph))
    if (not 0<a<1 or np.any(d<=0) or np.any(d>=1) or
        np.any(b<0) or not np.all(np.isfinite(d+b))):
        raise ValueError('physical Pauli/Bose occupations required')
    F=a*(1-d)*(1+b)
    R=(1-a)*d*b
    flux=F-R
    entropy=np.where(R>0,flux*np.log(np.maximum(F,1e-300)/np.maximum(R,1e-300)),
                     np.inf if np.any(F>0) else 0.)
    return {'mother_energy':e1,'daughter_E2':e2,'daughter_Ephi':eph,
            'forward':F,'reverse':R,'flux':flux,'entropy':entropy,
            'energy_defect':float(np.max(abs(e1-e2-eph)))}


def weak_moments(f1,f2,bphi,*,m1=.1,m2=.02,mphi=.001,y=.001,
                 pmax=12.,nr=12,nangle=16,g1=4):
    width=vacuum_width(m1,m2,mphi,y)
    if width<=0:return {'closed_channel':True,'gamma_rest':0.,
            'R_phi':0.,'R_chi2':0.,'R_chi1':0.,'Qtotal':0.,'entropy':0.}
    if not(pmax>0 and nr>=2 and nangle>=2 and g1>0):raise ValueError
    px,pw=leggauss(nr)
    ps=pmax*(px+1)/2
    ws=pw*pmax/2
    cc,wc=leggauss(nangle)
    total=heat1=heat2=heatphi=entropy=0.
    defect=0.
    for p,wp in zip(ps,ws):
        ev=event(float(p),cc,m1,m2,mphi,f1,f2,bphi)
        gamma=width*m1/ev['mother_energy']
        weight=g1/(2*math.pi**2)*wp*p*p*gamma/2
        total+=weight*np.dot(wc,ev['flux'])
        heat1-=weight*ev['mother_energy']*np.dot(wc,ev['flux'])
        heat2+=weight*np.dot(wc,ev['daughter_E2']*ev['flux'])
        heatphi+=weight*np.dot(wc,ev['daughter_Ephi']*ev['flux'])
        entropy+=weight*np.dot(wc,ev['entropy'])
        defect=max(defect,ev['energy_defect'])
    return {'closed_channel':False,'gamma_rest':width,
         'R_chi1':-total,'R_chi2':total,'R_phi':total,
         'Qchi1':heat1,'Qchi2':heat2,'Qphi':heatphi,
         'Qtotal':heat1+heat2+heatphi,
         'entropy':entropy,'max_event_energy_defect':defect,
         'synthetic_weak_only':True,'physical_relic_validated':False}


if __name__=='__main__':
 import json
 f1=lambda E:expit(-np.asarray(E)/.9)
 f2=lambda E:expit(-np.asarray(E))
 b=lambda E:1/np.expm1(np.asarray(E))
 print(json.dumps(weak_moments(f1,f2,b),indent=2))
