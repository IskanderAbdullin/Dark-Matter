"""Independent massive-boson on-shell thermal weak source with chemical potentials.

The two daughters share T0 but may have nonzero chemical potentials; the
parent has independent T1 and mu1. This is a prescribed-occupation reference,
NOT a strong kinetic solver, finite SM plasma or physical relic prediction.
"""
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.special import expit, spence

M1, M2, MPH, Y, G1 = .1, .02, 86e-9, .001, 4.


def kinematics(p):
    lam=(M1*M1-(M2+MPH)**2)*(M1*M1-(M2-MPH)**2)
    q=math.sqrt(lam)/(2*M1)
    estar=math.hypot(q,MPH)
    e1=math.hypot(p,M1)
    a=e1*estar/M1
    b=p*q/M1
    elo=(p*p*MPH*MPH+M1*M1*estar*estar)/(M1*(e1*estar+p*q))
    return e1,elo,a+b,b


def width():
    lam=(M1*M1-(M2+MPH)**2)*(M1*M1-(M2-MPH)**2)
    return Y*Y*((M1+M2)**2-MPH**2)*math.sqrt(lam)/(16*math.pi*M1**3)


def moments(e1, lo, hi, temp, mu2, muphi):
    """Exact integrals of (1,k)*fFD(e1-k;T,mu2)*bBE(k;T,muphi) dk."""
    assert temp>0 and muphi<MPH and 0<lo<=hi and e1-hi>=M2-1e-12
    a=math.exp((e1-mu2-muphi)/temp)
    def primitive(k):
        u=(k-muphi)/temp
        t=math.exp(u)
        logratio=math.log(-math.expm1(-u)) + u - float(np.logaddexp(u,math.log(a)))
        P=-spence(t) - u*math.log1p(t/a) - spence(1+t/a)
        return logratio,P
    lo0,lo1=primitive(lo)
    hi0,hi1=primitive(hi)
    i0=temp*(hi0-lo0)/(a+1)
    i1=muphi*i0+temp**2*(hi1-lo1)/(a+1)
    return i0,i1


def analytic(pmax, temp1, temp0, mu1, mu2, muphi):
    assert pmax>0 and temp1>0 and temp0>0
    def diff(p):
        e1,lo,hi,b=kinematics(p)
        f1=expit(-(e1-mu1)/temp1)
        affinity=-(e1-mu1)/temp1+(e1-mu2-muphi)/temp0
        I0,I1=moments(e1,lo,hi,temp0,mu2,muphi)
        if b == 0:
            inv=expit(-((e1-lo)-mu2)/temp0)/math.expm1((lo-muphi)/temp0)
            i0=inv;i1=lo*inv
        else:
            i0=I0/(2*b);i1=I1/(2*b)
        weight=G1*p*p/(2*math.pi**2)*width()*M1/e1
        j=weight*(1-f1)*math.expm1(affinity)
        return j*i0, j*i1, -j*e1*i0
    vals=[quad(lambda p: diff(p)[k],0,pmax,epsabs=1e-19,
               epsrel=2e-10,limit=240)[0] for k in (0,1,2)]
    return {'Rphi':vals[0], 'Qphi':vals[1], 'Q1':vals[2],
            'Q2':-vals[2]-vals[1]}


def gauss(pmax,nrad,nangle,temp1,temp0,mu1,mu2,muphi):
    x,w=leggauss(nrad); c,cw=leggauss(nangle)
    rate=qphi=q1=q2=h=0.
    for xx,ww in zip(x,w):
        p=pmax*(xx+1)/2;dp=ww*pmax/2
        e1,lo,hi,b=kinematics(p)
        if not b>0:raise ArithmeticError('nonpositive boosted interval')
        d=math.log(hi/lo)
        k=lo*np.exp((c+1)*d/2)
        dcos=cw*k*d/(2*b)
        W=G1*p*p*dp/(2*math.pi**2)*width()*M1/e1*dcos/2
        e2=e1-k
        f1=expit(-(e1-mu1)/temp1)
        f2=expit(-(e2-mu2)/temp0)
        bose=1/np.expm1((k-muphi)/temp0)
        affinity=-(e1-mu1)/temp1+(e1-mu2-muphi)/temp0
        j=W*(1-f1)*f2*bose*math.expm1(affinity)
        rate+=float(sum(j));qphi+=float(j@k)
        q2+=float(j@e2);q1-=float(sum(j)*e1)
        h+=float(sum(j)*affinity)
    return {'Rphi':rate,'Qphi':qphi,'Q1':q1,'Q2':q2,
            'energy_residual':q1+q2+qphi,'entropy_rate':h}


def audit():
    cases=[]
    for label,t1,mu1 in [('equilibrium',1.,-.15),
                         ('forward',1.,-.09),
                         ('reverse',1.,-.21),
                         ('temperature_split',.9,-.15)]:
        for pmax in (2.,8.,16.):
            ref=analytic(pmax,t1,1.,mu1,.05,-.2)
            numerical=gauss(pmax,72,48,t1,1.,mu1,.05,-.2)
            cases.append({'case':label,'pmax':pmax,
                          'analytic':ref,'gauss':numerical})
    return cases


if __name__=='__main__':
    import json
    print(json.dumps(audit(),indent=2,allow_nan=False))
