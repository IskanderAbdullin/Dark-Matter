"""Complete the vacuum Yukawa radiator over the exact hard real support."""
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre

from early_pole_hard_yukawa_photoproduction_reference import invariant_form
from finite_mass_qed_absolute_vacuum_window_reference import asym
from finite_mass_qed_signed_real_soft_recoil_reference import setup,angular_scaled


def acceptance(m,s,r):
    S=1-2*r
    p=np.sqrt(max(0.,S-(m+s)**2)*(S-(m-s)**2))/(2*np.sqrt(S))
    E=np.hypot(p,m)
    a=2*r/np.sqrt(S)
    return m*m+a*m*m/(E+p),m*m+a*(E+p)


def band_rest_angle(m,s,a,b,n, moment=None):
    z,w,p,E,_=setup(1.,m,s,n)
    norm=p*4*(E+m)
    def kernel(r):
        amp=angular_scaled(r,z,1.,m,s,np.sqrt)
        if moment=='scalar_loss':
            S=1-2*r
            q=np.sqrt((S-(m+s)**2)*(S-(m-s)**2))/(2*np.sqrt(S))
            charge=np.hypot(q,m)
            return (w@(amp*(m*m/(charge+q)+q*(1-z))))/(norm*np.sqrt(S))
        return (w@amp)/norm/(1. if moment=='photon' else r)
    def f(x):
        # Smoothly approach the exact upper threshold, including sqrt closure.
        return kernel(b-(b-a)*x*x)*2*(b-a)*x
    return quad(f,0,1,epsabs=8e-12,epsrel=5e-11,limit=140)[0]


def band_reversed_dalitz(m,s,a,b,n, moment=None):
    _,_,p,E,_=setup(1.,m,s,12)
    norm=p*4*(E+m)
    tmax=(1-s)**2
    nodes,weights=roots_legendre(n)
    stops=sorted({m*m,tmax,*acceptance(m,s,a),*acceptance(m,s,b)})
    def f(t):
        root=np.sqrt(t)
        Ep=(1+t-s*s)/(2*root)
        pstar=np.sqrt(((1-s)**2-t)*((1+s)**2-t))/(2*root)
        kg=(t-m*m)/(2*root)
        rp=Ep+pstar
        low=max(a,kg/rp)
        high=min(b,kg*rp)
        if high<=low:
            return 0.
        # r is positive in this branch. Integrate along log r, independently
        # of the parent-rest angular mapping at fixed photon energy.
        l,h=np.log(low),np.log(high)
        r=np.exp((l+h+(h-l)*nodes)/2)
        v=invariant_form(-r,-(t-m*m)/2,(1+t-s*s)/2-r,1.,m)*r
        if moment=='photon':
            v*=r
        if moment=='scalar_loss':
            v*=(t-m*m)/2
        return (h-l)*(weights@v)/(8*norm)
    integral=0.
    for l,h in zip(stops[:-1],stops[1:],strict=True):
        if h<=l:
            continue
        def mapped(x):
            return f(l+(h-l)*np.sin(x/2)**2)*np.sin(x)*(h-l)/2
        integral+=quad(mapped,0,np.pi,epsabs=6e-12,epsrel=5e-11,limit=200)[0]
    return integral


def compare(M,m,s):
    m,s=m/M,s/M
    kmax=(1-(m+s)**2)/2
    checks=[]
    for cap in (.001,.01,.05,.14,.71*kmax):
        if cap>=kmax:
            continue
        matched,ladder=asym(1.,m,s*s,cap)
        h1=band_rest_angle(m,s,cap,kmax,64)
        h2=band_reversed_dalitz(m,s,cap,kmax,56)
        h3=band_reversed_dalitz(m,s,cap,kmax,80)
        numerical=max(abs(h1-h3),abs(h2-h3))
        if numerical>2e-8:
            raise ArithmeticError(f'Independent hard-real phase spaces disagreed: {h1,h2,h3}.')
        checks.append(dict(separator=float(cap),unresolved_and_virtual=matched,hard_real=h3,full=matched+h3,
                           independent_real_defect=float(numerical),mass_regulator_last_defect=ladder[-1]['defect']))
    vals=[x['full'] for x in checks]
    spread=max(vals)-min(vals)
    inclusive,ladder=asym(1.,m,s*s,kmax)
    if spread>3e-8 or abs(vals[0]-inclusive)>4e-8 or ladder[-1]['defect']>5e-7:
        raise ArithmeticError(f'The full rate did not close under a changed artificial separator: {vals,inclusive}.')
    rad_energy=band_rest_angle(m,s,0,kmax,64,'photon')
    energy_ladder=[]
    for k in (1e-3,1e-5,1e-7):
        a=band_rest_angle(m,s,k,kmax,64,'photon')
        b=band_reversed_dalitz(m,s,k,kmax,64,'photon')
        if abs(a-b)>2e-8:
            raise ArithmeticError('Hard-real energy moment did not reproduce in reversed phase space.')
        energy_ladder.append(dict(minimum_k_over_M=k,finite_band_moment=b,independent_defect=abs(a-b),
                                  omitted_soft_energy=rad_energy-a))
    scalar_loss=band_rest_angle(m,s,0,kmax,64,'scalar_loss')
    scalar_band=band_rest_angle(m,s,1e-7,kmax,64,'scalar_loss')
    scalar_reverse=band_reversed_dalitz(m,s,1e-7,kmax,64,'scalar_loss')
    if abs(scalar_band-scalar_reverse)>2e-9 or not 0<scalar_loss<rad_energy:
        raise ArithmeticError('Scalar recoil energy moment did not match the reversed phase space.')
    return dict(endpoint_k_over_M=kmax,separator_checks=checks,
                stable_separator_spread=spread,absolute_inclusive_coefficient=inclusive,
                vacuum_rate_correction_at_g_005=inclusive*.005**2/(4*np.pi**2),
                real_radiated_energy_over_M_before_gauge=rad_energy,real_energy_ladder=energy_ladder,
                scalar_energy_shift_over_M_before_gauge=-scalar_loss,
                daughter_energy_shift_over_M_before_gauge=scalar_loss-rad_energy,
                scalar_energy_finite_band_defect=float(abs(scalar_reverse-scalar_band)),
                scalar_energy_omitted_below_1e7=float(scalar_loss-scalar_band))


def run():
    rows=[]
    for grid in ('fine','third'):
        path=Path(__file__).parent/f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz'
        with gzip.open(path) as f:d=json.load(f)
        M=d['x_end']
        m=M*d['mass_ratio']
        s=d['checkpoints'][-1]['scalar_pole_over_T']
        rows.append(dict(grid=grid,scalar_mass=s,**compare(M,m,s)))
    return dict(scope='O(g²) complete rest-frame vacuum decay with on-shell massive fermions and MS-bar Yukawa at mu=M; no thermal or cosmological operator.',rows=rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
