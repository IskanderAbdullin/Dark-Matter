"""Leading soft radiation and a measurement-level real/virtual counterterm."""
import gzip
import json
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.integrate import quad

from early_pole_decay_mode_clock_reference import Q, intervals
from early_pole_hard_yukawa_photoproduction_reference import invariant_form


def rest_geometry(M,m,ms):
    p=np.sqrt((M*M-(m+ms)**2)*(M*M-(m-ms)**2))/(2*M)
    E=np.hypot(m,p)
    chi=np.arcsinh(p/m)
    reference=2*(chi/np.tanh(chi)-1.)
    return p,E,chi,reference


def angular_rest(M,m,ms):
    p,E,chi,known=rest_geometry(M,m,ms)
    beta=p/E
    L,R=np.log((m/E)**2/(1+beta)),np.log1p(beta)
    def integrand(y):
        r=np.exp(y)
        z=(1-r)/beta
        return beta*.5*(1-z*z)/r
    return quad(integrand,L,R,epsabs=2e-13)[0],chi,known


def angular_lab(M,m,ms,q):
    p0,E0,_,ref=rest_geometry(M,m,ms)
    Es,l,h,*_=intervals(M,ms,m,q)
    Ep=(l+h)/2
    pp=np.sqrt((Ep-M)*(Ep+M))
    ga=Ep/M
    be=pp/Ep
    az=(Es/ga-np.hypot(p0,ms))/(be*p0)
    if abs(az)>1.:
        raise ValueError('Scalar mode did not intersect the vacuum aperture.')
    E=ga*(E0-be*p0*az)
    k=ga*(-p0*az+be*E0)
    d=M*E0
    # Cross term after azimuth integration. For z=cos(theta),
    # E²-kz²-k_perp²=m² gives a positive, subtracted discriminant.
    def fn(t):
        z=1-t
        pd=M*M/(Ep+pp)+pp*t
        radical=np.hypot(E*z-k,m*np.sqrt(t*(2-t)))
        return d/(pd*radical)
    x=np.r_[0.,np.geomspace(1e-8,.02,16),2.]
    total=sum(quad(fn,a,b,epsabs=2e-13,epsrel=2e-13)[0] for a,b in zip(x[:-1],x[1:],strict=True))-2
    return total,ref,float(pp)


def exact_soft_limit(M,m,ms):
    p0,E0,_,_=rest_geometry(M,m,ms)
    angles=np.r_[-.87,-.4,.2,.65,.91,.977]
    expected=(p0/E0)**2*(1-angles*angles)/(1-(p0/E0)*angles)**2
    out=[]
    for r in (1e-3,1e-5,1e-7,1e-9):
        k=M*r
        s=M*(M-2*k)
        mass=np.sqrt(s)
        p=np.sqrt((s-(m+ms)**2)*(s-(m-ms)**2))/(2*mass)
        ep=np.hypot(p,m)
        boost=-k/(M-k)
        ga=(M-k)/mass
        out_e=ga*(ep+boost*p*angles)
        out_z=ga*(p*angles+boost*ep)
        # Exact emitted photon signed invariant; not a soft approximation.
        exact=invariant_form(-M*k,-k*(out_e-out_z),M*out_e,M,m)
        scaled=exact*k*k/(4*M*(E0+m))
        error=float(max(abs(scaled/expected-1)))
        out.append(dict(rest_photon_energy_over_parent_mass=r,maximum_relative_soft_residual=error))
    if out[-1]['maximum_relative_soft_residual']>1e-6:
        raise ValueError('Coherent real matrix element failed the leading soft limit.')
    return out


def cancellation(temperature=.9):
    c=.7
    upper=.02
    mu0=1e-12
    def smooth(omega):
        B=1/np.expm1(omega/temperature)
        return (np.expm1(c*omega)+4*B*np.sinh(c*omega/2)**2)/omega
    rate=[]
    target=quad(smooth,0.,upper,epsabs=1e-13)[0]
    for cut in (1e-3,1e-5,1e-8,mu0):
        # Log-coordinate quadrature integrates the separate diverging pieces.
        lims=np.log([cut,upper])
        separate=quad(lambda y:1+2/np.expm1(np.exp(y)/temperature),*lims,
                      epsabs=1e-12,epsrel=1e-13)[0]
        stable=quad(smooth,cut,upper,epsabs=1e-13)[0]
        cusp=quad(lambda y:c*np.exp(y)*(1+2/np.expm1(np.exp(y)/temperature)),
                  *lims,epsabs=1e-13)[0]
        rate.append(dict(lower_cut=cut,real_plus_constant_measurement=separate,
                         leading_virtual_constant_measurement=-separate,
                         smooth_combined=stable,smooth_error_to_zero_cut=abs(target-stable),
                         nonsmooth_cusp_combined=cusp))
    with mp.workdps(80):
        T=mp.mpf(temperature)
        C=mp.mpf(c)
        def uncombined(w):
            N=1/mp.expm1(w/T)
            return ((1+N)*mp.exp(C*w)+N*mp.exp(-C*w)-(1+2*N))/w
        verify=float(mp.quad(uncombined,[mp.mpf(mu0),mp.mpf('.00001'),mp.mpf(upper)]))
    discrepancy=abs(verify-rate[-1]['smooth_combined'])
    cusp_slope=(rate[-1]['nonsmooth_cusp_combined']-rate[-2]['nonsmooth_cusp_combined'])/np.log(rate[-2]['lower_cut']/mu0)
    if discrepancy>2e-14 or abs(cusp_slope-2*c*temperature)>1e-8:
        raise ValueError('Thermal real/virtual cancellation or nonsmooth control failed.')
    return dict(smooth_zero_cut_result=target,high_precision_combination_defect=discrepancy,
                cusp_log_coefficient=cusp_slope,rows=rate)


def run():
    rows=[]
    for g in ('fine','third'):
        name=Path(__file__).parent/f'results/early-pole-spectrum-{g}-m06-2026-11-26.json.gz'
        with gzip.open(name) as f:d=json.load(f)
        M=d['x_end'];m=M*d['mass_ratio'];ms=d['checkpoints'][-1]['scalar_pole_over_T']
        a,c,e=angular_rest(M,m,ms)
        b,_,p=angular_lab(M,m,ms,Q)
        if max(abs(a-e),abs(b-e))>1e-8:
            raise ValueError('Leading soft angular kernel is not invariant across two tested frames.')
        rows.append(dict(grid=g,relative_rapidity=c,leading_radiator=e,independent_rest_angular=a,
                         independent_plasma_angular=b,plasma_parent_momentum=p,
                         matrix_element_soft_sequence=exact_soft_limit(M,m,ms)))
    return dict(scope='Leading-soft eikonal subtraction, artificial measurement weights, not matched physical finite-temperature radiation.',
                angular=rows,subtraction=cancellation())


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
