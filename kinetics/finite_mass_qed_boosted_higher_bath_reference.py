"""Cold photon-only current: boost-clock and higher even recoil checks."""
import gzip
import json
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre

from finite_mass_qed_analytic_low_t_identity_reference import exact_real


def real_primitive(m,s,current):
    m,s=map(mp.mpf,(m,s))
    u=(1+m*m-s*s)/2
    p=mp.sqrt((1-(m+s)**2)*(1-(m-s)**2))/2
    a=m if current=='scalar' else 0
    def f(r):
        S=1-2*r
        hi=(S+m*m-s*s+mp.sqrt((S-(m+s)**2)*(S-(m-s)**2)))/(2*S)
        lo=m*m/(S*hi)
        def primitive(w):
            return -4*(u+a)*(w-m*m/w-2*u*mp.log(w))+8*r*(u+a)*(w-mp.log(w))+r*r*(2*w*w-8*w+4*mp.log(w))
        return (primitive(hi)-primitive(lo))/(8*p*(u+a))
    return f


def doppler_moment(beta,n):
    # Integration over cos(theta) of [gamma(1+beta*cos(theta))]^-n.
    if beta==0:
        return 1.
    ell=np.arctanh(beta)
    gamma=np.cosh(ell)
    return np.sinh((n-1)*ell)/(gamma*beta*(n-1))


def average_occupation(t,beta):
    # t is measured at the minimum Doppler exponent in the angular domain.
    if beta==0:
        return np.exp(-t)/(-np.expm1(-t))
    ratio=(1+beta)/(1-beta)
    return (np.log1p(-np.exp(-ratio*t))-np.log(-np.expm1(-t)))/(t*(ratio-1))


def geometry_check(m,s,beta):
    u=(1+m*m-s*s)/2
    p=np.sqrt((u-m)*(u+m))
    daughter_exact=(np.log(u+p)-np.log(m*m/(u+p)))/p
    gamma=1/np.sqrt((1-beta)*(1+beta))
    x,w=roots_legendre(92)
    if beta:
        l,r=np.log1p(-beta),np.log1p(beta)
        y=(l+r+(r-l)*x)/2
        cu=(np.exp(y)-1)/beta
        wp=(r-l)*np.exp(y)*w/(4*beta)
    else:
        cu,wp=x,w/2
    D=gamma*(1+beta*cu)
    Jy=(np.log(m*m/(u+p)),np.log(u+p))
    ld=(Jy[0]+Jy[1]+(Jy[1]-Jy[0])*x)/2
    cd=(u-np.exp(ld))/p
    wd=(Jy[1]-Jy[0])*np.exp(ld)*w/(4*p)
    nk=np.stack([np.sqrt(1-cu*cu),0*cu,cu],axis=1)
    ex=np.stack([cu,0*cu,-np.sqrt(1-cu*cu)],axis=1)
    az=2*np.pi*np.arange(12)/12
    nb=nk[:,None,None,:]*cd[None,:,None,None]+np.sqrt(1-cd*cd)[None,:,None,None]*(
        ex[:,None,None,:]*np.cos(az)[None,None,:,None]+np.eye(3)[1]*np.sin(az)[None,None,:,None])
    dot=np.sum(nk[:,None,None,:]*nb,axis=-1)
    gap=m*m/(u+p)+p*(1-dot)
    dip=wp/D**2
    parent=2*(dip@(gap.mean(axis=2)@wd))
    daughter=2*(dip@((1/gap).mean(axis=2)@wd))
    spin_def=max(abs(parent/(2*u)-1),abs(daughter/daughter_exact-1))
    # Independent angular averaging and analytic Bose occupation moment.
    mn=gamma*(1-beta)
    moment_def=[]
    for n in (2,4,6):
        Y=quad(lambda t:t**(n-1)*average_occupation(t,beta),0,67,
               epsabs=3e-11,epsrel=4e-11)[0]/mn**n
        rhs=float(mp.gamma(n)*mp.zeta(n))*doppler_moment(beta,n)
        moment_def.append(max(abs(Y/rhs-1),abs((wp@(1/D**n))/doppler_moment(beta,n)-1)))
    # T^2 comparison in the bath clock at fixed spatial momentum, not
    # at fixed velocity. The normalized proper rate is 1-2*h.
    with mp.workdps(45):
        pl=mp.mpf(gamma)*mp.mpf(beta)
        g0=mp.sqrt(1+pl**2)
        def rate_for_h(h):
            mt=mp.sqrt(1+4*h)
            return (1-2*h)*g0*mt/mp.sqrt(pl**2+mt**2)
        via_clock=float(mp.diff(rate_for_h,0))
    if max(spin_def,max(moment_def),abs(via_clock+2/gamma**2))>2e-8:
        raise ArithmeticError('Boosted angular projection or clock replay failed.')
    return dict(beta=beta,gamma=gamma,spin_trace_relative_defect=float(spin_def),
                doppler_average_inverse_2=doppler_moment(beta,2),
                doppler_average_inverse_4=doppler_moment(beta,4),
                doppler_average_inverse_6=doppler_moment(beta,6),
                bose_and_angle_relative_defects=[float(z) for z in moment_def],
                proper_clock_T2_coefficient=-2.,
                bath_clock_fixed_momentum_T2_coefficient=via_clock)


def higher_controls(m,s,current,betas):
    with mp.workdps(58):
        f=real_primitive(m,s,current)
        soft=f(0)
        endpoint=(1-(m+s)**2)/2
        if max(abs(exact_real(m,s,current,r)-float(f(r))) for r in (-.15,.015,.65*endpoint))>3e-10:
            raise ArithmeticError('High-precision amplitude and public primitive did not match.')
        coeff={n:2*mp.diff(f,0,n)/mp.factorial(n) for n in (2,4,6,8)}
        radii=[mp.mpf(endpoint)/8,mp.mpf(endpoint)/12]
        roots=[mp.exp(2j*mp.pi*j/48) for j in range(48)]
        check={}
        for n in (4,6,8):
            c=[2*mp.re(mp.fsum((f(rad*z)-soft)*z**(-n) for z in roots)/48/rad**n) for rad in radii]
            check[n]=max(float(abs(v-coeff[n])) for v in c)
        if max(check.values())>1e-25:
            raise ArithmeticError('Independent complex contour did not resolve the high even coefficients.')
        receipts=[]
        for beta in betas:
            Dmin=np.sqrt((1-beta)/(1+beta))
            gamma=1/np.sqrt((1-beta)*(1+beta))
            nints={n:float(coeff[n]*mp.gamma(n)*mp.zeta(n))*doppler_moment(beta,n)*Dmin**n
                   for n in (4,6,8)}
            rows=[]
            for sigma in (.0004,.0008,.0016):
                tau=Dmin*sigma
                kmax=mp.mpf(endpoint)
                SM=mp.mpf(sigma)
                def integrand(t):
                    if not t:
                        return 0.
                    r=mp.mpf(t)*SM
                    residual=f(-r)+(f(r) if r<kmax else 0)-2*soft-coeff[2]*r*r
                    return float(residual/(SM**4*mp.mpf(t)))*average_occupation(t,beta)
                def integrate(top):
                    points=[1.,5.,10.,25.]
                    if endpoint/sigma<top:
                        points.append(endpoint/sigma)
                        points.sort()
                    return quad(integrand,0,top,points=points,epsabs=4e-12,epsrel=4e-11)[0]
                value=integrate(61.)
                plan=nints[4]+sigma*sigma*nints[6]+sigma**4*nints[8]
                discrepancy=abs(value-plan)
                tail=abs(value-integrate(73.))
                if tail>3e-9 or discrepancy>8e-8*max(1.,abs(value)):
                    raise ArithmeticError(f'The full Bose residual failed its power/endpoint replay: {value,plan,tail}.')
                rows.append(dict(max_directional_T_over_M=sigma,bath_temperature_over_M=tau,
                                 directly_integrated_higher_recoil_over_sigma4=value,
                                 expansion_through_T8_over_sigma4=plan,T8_replay_defect=discrepancy,
                                 independent_photon_tail_defect=tail))
            receipts.append(dict(beta=beta,gamma=gamma,leading_T4_before_bose_moment=float(coeff[4]),
                                 boosted_T4_in_fixed_bath_temperature=float(coeff[4])*doppler_moment(beta,4),
                                 rows=rows))
        return dict(m_over_M=m,scalar_mass_over_M=s,current=current,
                    signed_even_recoil_coefficients={n:float(c) for n,c in coeff.items()},
                    independent_contour_defects=check,boosted_bose_checks=receipts)


def run():
    outputs=[]
    for name in ('fine','third'):
        path=Path(__file__).parent/f'results/early-pole-spectrum-{name}-m06-2026-11-26.json.gz'
        with gzip.open(path) as f:d=json.load(f)
        m=d['mass_ratio'];s=d['checkpoints'][-1]['scalar_pole_over_T']/d['x_end']
        outputs.append(dict(mass_fixture=name,**higher_controls(m,s,'scalar',(0.,.8,.985))))
    outputs.extend(higher_controls(m,s,current,(0.,.8))
                   for m,s,current in ((.2,.65,'scalar'),(.2,.65,'left'),(.5,.12,'scalar'),(.5,.12,'left')))
    return dict(scope='Cold bare photon-only isotropic bath in its rest frame: moving equal charges, no plasma dispersion, charged occupations or hot-archive claim.',
                controls=outputs,
                independent_spin_and_clock=[geometry_check(.06,.0174,b) for b in (0.,.4,.8,.95,.985,.998)])


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
