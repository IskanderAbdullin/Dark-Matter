"""Formal photon-only comparison across the exact real-emission endpoint."""
import gzip
import json
from math import fsum
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar
from scipy.special import roots_legendre

from finite_mass_qed_analytic_low_t_identity_reference import exact_real
from finite_mass_qed_boosted_higher_bath_reference import average_occupation, real_primitive


def integrated_independent(m,s,current,r,n=64):
    S=1-2*r
    M2=m*m
    sr=np.sqrt((S-(m+s)**2)*(S-(m-s)**2))
    b=(S+M2-s*s+sr)/(2*S)
    a=M2/(S*b)
    u=(1+M2-s*s)/2
    p=np.sqrt((1-(m+s)**2)*(1-(m-s)**2))/2
    spin=m if current=='scalar' else 0.
    x,w=roots_legendre(n)
    l,h=np.log(a),np.log(b)
    eta=np.exp((l+h+(h-l)*x)/2)
    d=u+r*(eta-1)
    h0=1+M2/eta**2-2*d/eta
    v=(-1+1/eta)/2
    x0=eta+M2/eta-2*d
    kernel=-4*(d+spin)*h0-8*r*v*x0+16*r*r*v*v*eta
    return (h-l)/2*(w@(kernel*eta))/(8*p*(u+spin))


def fixture(m,s,current):
    cap=(1-(m+s)**2)/2
    with mp.workdps(50):
        f=real_primitive(m,s,current)
        origin=float(f(0))
        c2,c4,c6=(float(2*mp.diff(f,0,j)/mp.factorial(j)) for j in (2,4,6))
    samples=[-.001,-.12,-.8,-6.,-30.,-240.,-.0002,.0006,.4*cap,.98*cap]
    defects=[abs(integrated_independent(m,s,current,r)/exact_real(m,s,current,r)-1) for r in samples]
    if max(defects)>5e-10:
        raise ArithmeticError('Wide real angular calculation did not reproduce the exact primitive.')
    def compute(tau,beta,direct=False):
        loD=np.sqrt((1-beta)/(1+beta))
        scale=tau/loD
        Tcap=cap/scale
        qlimits=sorted({0.,.0005/scale,Tcap,1.,4.,10.,24.,48.,74.})
        qlimits=[x for x in qlimits if x<=74.]
        def f(t):
            r=t*scale
            if r<.0005:
                w=r*r
                G=w*(c2+w*(c4+w*c6))
            else:
                kernel=integrated_independent if direct else exact_real
                G=kernel(m,s,current,-r)+(kernel(m,s,current,r) if r<cap else 0.)-2*origin
            return G*average_occupation(t,beta)/t
        if direct:
            x,w=roots_legendre(34)
        ss=[]
        for a,b in zip(qlimits[:-1],qlimits[1:],strict=True):
            if direct:
                phi=(x+1)*np.pi/2
                u=a+(b-a)*np.sin(phi/2)**2
                vals=np.array([f(z) for z in u])
                result=((w*np.sin(phi))@vals)*(b-a)*np.pi/4
            else:
                result=quad(f,a,b,epsabs=3e-11,epsrel=2e-10,limit=170)[0]
            ss.append((a>=Tcap,result))
        J=np.pi**2/6*tau*tau
        below=fsum(v for k,v in ss if not k)/J
        above=fsum(v for k,v in ss if k)/J
        # On-shell external leading projection in the matched convention.
        external=-2-c2
        proper=below+above+external
        return np.array([proper,below,above,external,proper+2*beta*beta])
    return cap,max(defects),compute


def run_cases(m,s,current,betas,find_roots):
    cap,angular,calc=fixture(m,s,current)
    results=[]
    for beta in betas:
        checks=[]
        for T in (.002,.01,.04,.1,.2,.4,.8,1.6,3.2,6.4,12.8,20.):
            a=calc(T,beta)
            if T in (.002,.1,.8,6.4,20.):
                d=max(abs(a-calc(T,beta,True)))
                if d>3e-7:
                    raise ArithmeticError(f'Independent all-support thermal replay failed: {d}.')
            else:
                d=None
            checks.append(dict(T_over_M=T,proper_clock_coefficient=a[0],
                               signed_real_below_emission_endpoint=a[1],signed_real_above_endpoint=a[2],
                               external_mass_and_spin=a[3],bath_clock_coefficient_at_fixed_p=a[4],
                               independent_integral_defect=d))
        r=[]
        for a,b in zip(checks[:-1],checks[1:],strict=True):
            if a['proper_clock_coefficient']*b['proper_clock_coefficient']<0:
                r.append(brentq(lambda x:calc(x,beta)[0],a['T_over_M'],b['T_over_M'],xtol=2e-9))
        extreme=None
        if find_roots:
            fit=minimize_scalar(lambda x:calc(np.exp(x),beta)[0],bounds=(np.log(.2),np.log(20.)),
                                method='bounded',options={'xatol':2e-8})
            extreme=dict(T_over_M=float(np.exp(fit.x)),coefficient=float(fit.fun))
        results.append(dict(beta=beta,controls=checks,bracketed_proper_clock_sign_changes=r,
                            sampled_proper_clock_minimum=extreme))
    return dict(m_over_M=m,s_over_M=s,current=current,emission_endpoint_k_over_M=cap,
                independent_wide_angular_relative_defect=angular,moving=results)


def run():
    outputs=[]
    for name in ('fine','third'):
        p=Path(__file__).parent/f'results/early-pole-spectrum-{name}-m06-2026-11-26.json.gz'
        with gzip.open(p) as f:d=json.load(f)
        m=d['mass_ratio'];s=d['checkpoints'][-1]['scalar_pole_over_T']/d['x_end']
        outputs.append(dict(mass_fixture=name,**run_cases(m,s,'scalar',(0.,.8),True)))
    outputs.extend(run_cases(.2,.65,k,(0.,),True) for k in ('scalar','left'))
    return dict(scope='Formal O(g²) bare photon-only comparison after IR matching. At large T this is not a controlled entire plasma or an archive state; omits charged bath, widths, and photon dispersion.',
                units='Divide relative-to-Born correction by [g²/(4 pi²)] (pi²/6) (T/M)^2.',
                rows=outputs)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
