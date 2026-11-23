"""Exact boosted daughter energy log-angle quadrature for χ1↔χ2+real φ.

A cold Bose daughter requires resolution of Ephi=a-b*cos(theta*),
which may be sharply peaked near cos(theta*)=+1. y=log(Ephi) gives
positive dc = Ephi/b dy. This is a purely numerical on-shell 1↔2
quadrature, not a physical Bose condensate or cosmological rate.
"""
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import expit
from kinetics.phi_bose_cell_entropy import kinematics,scalar_yukawa_width


def angular_nodes(p,m1,m2,mphi,n=12,scheme='log_ephi'):
    """Positive nodes, weights in CM rest-frame dcos(theta*)."""
    if not(p>=0 and m1>m2+mphi>=0 and n>=2):
        raise ValueError('open physical daughter kinematics')
    c,w=leggauss(n)
    if scheme=='uniform_cos' or p==0:
        return c,w
    if scheme!='log_ephi':raise ValueError('unknown angle map')
    E_min=kinematics(p,1.,m1,m2,mphi)[2]
    E_max=kinematics(p,-1.,m1,m2,mphi)[2]
    if not(E_min>0 and E_max>E_min):raise ValueError('nonpositive boosted Bose daughter energy')
    b=(E_max-E_min)/2
    # y=log(Ephi/Emin) fixes the logarithmic Rayleigh-Jeans 1/Ephi integral exactly.
    L=math.log(E_max/E_min)
    Ephi=E_min*np.exp(L*(c+1)/2)
    out_c=1-(Ephi-E_min)/b
    weights=w*L/2*Ephi/b
    return out_c,weights


def analytic_angular_moments(p,m1,m2,mphi,n=12,scheme='log_ephi'):
    c,w=angular_nodes(p,m1,m2,mphi,n=n,scheme=scheme)
    E=np.array([kinematics(p,float(x),m1,m2,mphi)[2] for x in c])
    e_min=kinematics(p,1.,m1,m2,mphi)[2]
    e_max=kinematics(p,-1.,m1,m2,mphi)[2]
    b=(e_max-e_min)/2
    return {'p':p,'n':n,'scheme':scheme,'measure':float(sum(w)),
            'RJ_1_over_E':float(np.dot(w,1/E)),
            'exact_measure':2.,'exact_RJ_1_over_E':math.log(e_max/e_min)/b}


def physical_events(*,m1=.1,m2=.02,mphi=86e-9,yukawa=.001,pmax=2.,nr=24,nangle=12,scheme='log_ephi',g1=4):
    if not (pmax>0 and nr>=2 and nangle>=2):raise ValueError
    x,wp=leggauss(nr)
    width=scalar_yukawa_width(m1,m2,mphi,yukawa)
    events=[]
    for v,w in zip(x,wp):
        p=pmax*(v+1)/2
        measure=g1/(2*math.pi**2)*(w*pmax/2)*p*p*width*m1/math.hypot(p,m1)
        cc,ww=angular_nodes(float(p),m1,m2,mphi,nangle,scheme)
        events.extend((float(p),float(c),float(measure*a/2)) for c,a in zip(cc,ww))
    return events,width


def weak_rate(shape='thermal',*,pmax=2.,nr=24,nangle=12,scheme='log_ephi'):
    """Independent continuous-occupation weak Rphi and entropy, same on-shell tree."""
    m1=.1;m2=.02;mphi=86e-9
    events,width=physical_events(m1=m1,m2=m2,mphi=mphi,pmax=pmax,nr=nr,
                                  nangle=nangle,scheme=scheme)
    Rphi=H=0.
    for p,c,W in events:
        E1,E2,Ephi=kinematics(p,c,m1,m2,mphi)
        eta=Ephi
        if shape=='broad':eta+=.30*math.exp(-((Ephi-.10)/.06)**2)
        elif shape=='sharp':eta+=.40*math.exp(-((Ephi-.04)/.02)**2)
        elif shape!='thermal':raise ValueError(shape)
        f1=expit(-E1/.9);f2=expit(-E2);b=1/math.expm1(eta)
        F=f1*(1-f2)*(1+b)
        I=(1-f1)*f2*b
        J=F-I
        Rphi+=W*J
        H+=W*J*(-E1/.9+E2+eta)
    return {'Rphi':Rphi,'H':H,'width':width,'pmax':pmax,'nr':nr,
            'nangle':nangle,'angle_scheme':scheme,'shape':shape}


if __name__=='__main__':
    import json
    for s in ('thermal','broad','sharp'):
        for n in (8,12,24,48):
            print(json.dumps(weak_rate(s,nangle=n),indent=2))
