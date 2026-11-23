"""Independent exact-angle Bose/Yukawa thermal weak-number reference.

For FD chi1 at T1 != T2=Tphi=T0 and mu=0 the on-shell
forward/reverse affinity is E1*(1/T0-1/T1), independent of daughter
angle. An exact analytic FD(E1-Ephi)*BE(Ephi) energy primitive
reduces the physical source to one parent-momentum integral.
The separately written positive event quadrature checks number,
energy and entropy. This is not an evolving physical relic solver.
"""
import json
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.special import expit

M1, M2, MPH, YUKAWA, G1, T1, T0 = .1, .02, 86e-9, .001, 4., .9, 1.


def daughter_bounds(p):
    lam=(M1*M1-(M2+MPH)**2)*(M1*M1-(M2-MPH)**2)
    q=math.sqrt(lam)/(2*M1)
    estar=math.hypot(q,MPH)
    E1=math.hypot(p,M1)
    a=E1*estar/M1
    b=p*q/M1
    elo=(p*p*MPH*MPH+M1*M1*estar*estar)/(M1*(E1*estar+p*q))
    ehi=a+b
    return E1,elo,ehi,b


def width():
    lam=(M1*M1-(M2+MPH)**2)*(M1*M1-(M2-MPH)**2)
    return YUKAWA**2*((M1+M2)**2-MPH**2)*math.sqrt(lam)/(16*math.pi*M1**3)


def _log_expm1(v):
    return v+math.log(-math.expm1(-v))


def scalar_thermal_number_integral(E1, elo, ehi):
    """Exact positive integral FD(E1-k;T0)*BE(k;T0) dk."""
    if not(0<elo<=ehi<=E1-M2+1e-12):raise ValueError('on-shell positive daughters')
    def L(k):
        z=k/T0
        return _log_expm1(z)-float(np.logaddexp(z,E1/T0))
    return T0*expit(-E1/T0)*(L(ehi)-L(elo))


def parent_differential(p):
    E1,elo,ehi,b=daughter_bounds(p)
    prefactor=G1/(2*math.pi**2)*p*p*width()*M1/E1
    affinity=E1*(1/T0-1/T1)
    inverse=(1-expit(-E1/T1))*scalar_thermal_number_integral(E1,elo,ehi)
    if p==0:
        f2=expit(-(E1-elo)/T0)
        bos=1/math.expm1(elo/T0)
        angular_inverse=f2*bos
    else:
        angular_inverse=inverse/(1-expit(-E1/T1))/(2*b)
    return prefactor*(1-expit(-E1/T1))*math.expm1(affinity)*angular_inverse


def analytic_number_rate(pmax):
    if not pmax>0:raise ValueError('positive parent cutoff')
    result,error=quad(parent_differential,0,pmax,epsabs=1e-17,epsrel=2e-11,limit=240)
    return {'Rphi':result,'quad_abserr':error,'pmax':pmax}


def gauss_weak(pmax,nrad,nangle):
    xs,ws=leggauss(nrad)
    cx,cw=leggauss(nangle)
    rate=Qphi=H=Q1=Q2=0.
    min_weight=math.inf
    min_energy=math.inf
    for x,w in zip(xs,ws):
        p=pmax*(x+1)/2
        dp=w*pmax/2
        E1,elo,ehi,b=daughter_bounds(p)
        if not b>0:raise ArithmeticError('positive boosted interval')
        LL=math.log(ehi/elo)
        ephi=elo*np.exp((cx+1)*LL/2)
        dcos=cw*(LL/2)*ephi/b
        W=G1/(2*math.pi**2)*dp*p*p*width()*M1/E1*dcos/2
        e2=E1-ephi
        f1=expit(-E1/T1)
        f2=expit(-e2/T0)
        bphi=1/np.expm1(ephi/T0)
        min_weight=min(min_weight,float(min(W)))
        min_energy=min(min_energy,float(min(ephi)))
        affinity=E1*(1/T0-1/T1)
        flux=(1-f1)*f2*bphi*math.expm1(affinity)
        C=W*flux
        rate+=float(sum(C))
        Qphi+=float(C@ephi)
        Q2+=float(C@e2)
        Q1-=float(sum(C)*E1)
        H+=float(sum(C)*affinity)
    return {'Rphi':rate,'Qphi':Qphi,'Q2':Q2,'Q1':Q1,
            'total_energy_defect':Q1+Q2+Qphi,'H':H,
            'mean_Ephi_signed':Qphi/rate,'min_positive_event_weight':min_weight,
            'min_sampled_Ephi_GeV':min_energy,'pmax':pmax,'nrad':nrad,'nangle':nangle}


def audit():
    independent=[analytic_number_rate(p) for p in (2.,8.,16.)]
    cases=[gauss_weak(p,nr,na) for p in (2.,8.,16.)
           for nr,na in ((24,12),(48,24),(96,48))]
    reference={v['pmax']:v['Rphi'] for v in independent}
    for v in cases:
        v['relative_number_error_vs_analytic_1D']=v['Rphi']/reference[v['pmax']]-1
    return {'schema_version':1,'analytic_number_rate':independent,
            'on_shell_positive_quadratures':cases,
            'mass_and_temperatures_GeV':{'m1':M1,'m2':M2,'mphi':MPH,
                                         'T1':T1,'T2':T0,'Tphi':T0},
            'scope':'Exact independent thermal number source with positive on-shell event weights, physical weak energy and entropy, not a physical strong time integrator.',
            'class_ready':False}


if __name__=='__main__':
    print(json.dumps(audit(),indent=2,allow_nan=False))
