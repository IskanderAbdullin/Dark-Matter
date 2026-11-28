"""Resolve the angular PV as a bare charged loop approaches coincidence."""
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad

from early_pole_decay_mode_clock_reference import intervals,Q
from finite_mass_qed_thermal_charged_pv_geometry_reference import chart

D=np.longdouble


def hilbert(c,a,u):
    """PV integral 1/[(mu-u) sqrt((mu-c)^2+a^2)] on [-1,1]."""
    h=(u-c)/a
    root=np.sqrt(1+h*h)
    y=np.arcsinh(h)
    pp,mm=np.exp(y),np.exp(-y)
    def f(x):
        v=np.exp(np.arcsinh((x-c)/a))
        return np.log(abs(v-pp))-np.log(v+mm)
    return (f(1)-f(-1))/(a*root)


def solve_axes(M,m,ms,e,which):
    M,m,ms,e,q=map(D,(M,m,ms,e,Q))
    E=np.hypot(q,ms)
    P=np.sqrt((e-M)*(e+M))
    z=(e*E-(M*M+ms*ms-m*m)/2)/q
    if which=='parent':
        ma,mb,et,zz,sg=M,m,e,z,D(1)
    else:
        ma,mb,et,zz,sg=m,M,e-E,z-q,D(-1)
    p=np.sqrt((et-ma)*(et+ma))
    c=zz/p
    d=(M*M+m*m-ms*ms)/2
    C=(M+m)**2/(d+M*m)
    Cq=4-2*(M+m)*ma/(d+M*m)
    A0=ma*np.sqrt(1-c*c)/(p*p)
    H0=(E/q-c*et/p)/p
    target=2*d*np.arcsinh(H0/A0)/(sg*q*p**3*np.hypot(A0,H0))

    def full(delta):
        delta=D(delta)
        er=et+delta
        r=np.sqrt((er-ma)*(er+ma))
        g=2*sg*q*r
        # Difference from the exact external direction, not large cancellations.
        u=c+delta*(E/q-c*(er+et)/(r+p))/r
        rad=2*ma*ma*np.sinh((np.arcsinh(p/ma)-np.arcsinh(r/ma))/2)**2
        pp=p*r
        a=np.sqrt(2*rad*(1-c*c)/pp)
        if not -1<u<1 or a<=0:
            raise ValueError('Local charged overlap not inside the PV chart.')
        def residual_t(x):
            du=D(x)-c
            Y=pp*pp*(du*du+a*a)
            N=2*rad*pp*c*du-rad*rad
            X=Y-N
            yl,xl=np.sqrt(Y),np.sqrt(X)
            return -N/(2*xl*yl*(xl+yl))
        stage=sorted({-1.,1.,float(u),float(c),*[float(c+j*a) for j in (-64,-8,-1,1,8,64)]})
        stage=[z for z in stage if -1<=z<=1]
        f0=residual_t(float(u))
        # Bounded remainder after extracting the full local singular kernel.
        def div(t):
            v=D(t)
            if abs(v-u)<D(1e-13):
                p=D(1e-11)
                return float((residual_t(u+p)-residual_t(u-p))/(2*p))
            return float((residual_t(v)-f0)/(v-u))
        pv_rem=quad(div,-1.,1.,points=stage[1:-1],epsabs=2e-10,epsrel=2e-8,limit=300)[0]
        log=np.log((1-u)/(1+u))
        pv_ph=-hilbert(c,a,u)/(2*pp)+pv_rem+f0*log
        photons=-(np.arcsinh((1-c)/a)-np.arcsinh((-1-c)/a))/(2*pp)
        photons+=quad(lambda t:float(residual_t(t)),-1.,1.,points=stage[1:-1],
                      epsabs=2e-10,epsrel=1e-8,limit=300)[0]
        whole=(C*log+4*d*pv_ph)/g+Cq*photons/2
        return float(whole),float((u-c)/delta),float(a/abs(delta))

    return full,float(target),float(Cq/(2*p*p))


def run():
    out=[]
    for grid in ('fine','third'):
        src=Path(__file__).parent/f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz'
        with gzip.open(src) as f:d=json.load(f)
        M=d['x_end'];m=M*d['mass_ratio'];s=d['checkpoints'][-1]['scalar_pole_over_T']
        _,lo,hi,*_=intervals(M,s,m,Q)
        central=(lo+hi)/2
        for cut in ('parent','daughter'):
            style,target,log_coefficient=solve_axes(M,m,s,central,cut)
            raw=chart(M,m,s,central,cut)
            rows=[]
            for r in (1e-4,1e-5,1e-6,1e-7,1e-8):
                z=[]
                for sg in (-1,1):
                    v,b,a=style(sg*r)
                    cv=None
                    if r>=1e-5:
                        old=raw(central+sg*r)['principal_value_coefficient']
                        cv=abs(v/old-1)
                        if cv>2e-7:
                            raise ValueError('Analytic singular subtraction did not replay the ordinary PV at resolved shifts.')
                    z.append(dict(sign=sg,delta_times_principal_value=sg*r*v,shift_slope=b,width_coefficient=a,
                                  direct_angular_relative_defect=cv))
                odd=(z[0]['delta_times_principal_value']+z[1]['delta_times_principal_value'])/2
                even=(z[1]['delta_times_principal_value']-z[0]['delta_times_principal_value'])/(2*r)
                rows.append(dict(abs_energy_shift=r,sides=z,odd_coefficient=odd,even_part=even,even_after_log_subtraction=float(even-log_coefficient*np.log(r)),
                                 odd_limit_relative_defect=abs(odd/target-1)))
            if rows[-1]['odd_limit_relative_defect']>2e-5:
                raise ValueError('Numerical bare overlap did not reproduce the leading odd singularity.')
            observed=(rows[-1]['even_part']-rows[-2]['even_part'])/np.log(.1)
            logerr=abs(observed/log_coefficient-1)
            if logerr>.03:
                raise ValueError('Resolved even part did not recover its integrable logarithm.')
            out.append(dict(grid=grid,cut=cut,analytic_odd_limit=target,
                            expected_even_log_coefficient=log_coefficient,measured_even_log_coefficient=float(observed),
                            relative_even_log_slope_defect=float(logerr), rows=rows))
    return dict(scope='Local bare geometrical PV overlap asymptote only, no physical cut weights or residues.',
                rows=out)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
