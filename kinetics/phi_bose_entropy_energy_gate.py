"""Conditional ideal-gas entropy/energy bound, no collision closure assumed.
Run: python kinetics/phi_bose_entropy_energy_gate.py (requires scipy).
Units eV, hbar=c=kB=1.
"""
import json
from math import pi, sqrt, exp, expm1, log1p, log
from scipy.integrate import quad
from scipy.optimize import brentq

SPECIES = ((86., 1, False), (20e6, 4, True), (100e6, 4, True))
A = 2*pi*pi/45
NU = 7*pi*pi/120
HBARC = 1.973269804e-5  # eV cm
S0 = 2891.2  # cm^-3
RHO_CRIT = 1.05375e4  # eV cm^-3, rho_crit/h^2


def moments(t, species=SPECIES, log_pressure=False):
    """Return physical rho,P,s at zero chemical potentials."""
    rho = pressure = 0.
    for m,g,fd in species:
        z = m/t
        def f(y):
            e = sqrt(y*y+z*z)
            q = exp(-e)
            return q/(1+q) if fd else q/(-expm1(-e))
        # Upper cutoff 100 in p/T; exp(-100) tail for relativistic gas.
        coeff = g*t**4/(2*pi*pi)
        rho += coeff*quad(lambda y:y*y*sqrt(y*y+z*z)*f(y),0,100,epsabs=1e-30,epsrel=1e-9)[0]
        if log_pressure:
            def p_log(y):
                q=exp(-sqrt(y*y+z*z))
                return y*y*(log1p(q) if fd else -log1p(-q))
            pressure += coeff*quad(p_log,0,100,epsabs=1e-30,epsrel=1e-9)[0]
        else:
            pressure += coeff*quad(lambda y:y**4/(3*sqrt(y*y+z*z))*f(y),0,100,epsabs=1e-30,epsrel=1e-9)[0]
    return rho,pressure,(rho+pressure)/t


def min_energy(s, guess):
    t=brentq(lambda t:moments(t)[2]/s-1,guess/1e5,guess*10,xtol=1e-12,rtol=1e-11)
    r,p,st=moments(t)
    p2=moments(t,log_pressure=True)[1]
    assert abs(p2/p-1)<1e-7
    assert abs((t*s-p)/r-1)<1e-7
    return t,r,p


def main():
    tdec=.6882487043e9
    gd=moments(tdec)[2]/(A*tdec**3)
    r=gd/61.75
    result={'g_dark_dec':gd,'s_dark_over_s_sm':r,'bbn':[],'entropy_targets':[]}
    # Tref defined by total SM entropy; equals neutrino T in ideal instantaneous decoupling after e+e- annihilation.
    for t_ref in (3e6,2e6,1e6,.5e6,.2e6,.1e6,.05e6):
        s=r*A*10.75*t_ref**3
        t,rho,p=min_energy(s,t_ref)
        t_trial=(s/A)**(1/3)
        pt=moments(t_trial,log_pressure=True)[1]
        lower=t_trial*s-pt
        result['bbn'].append({'Tref_MeV':t_ref/1e6,'Tmin_over_Tref':t/t_ref,'rho_over_one_nu':rho/(NU*t_ref**4),'trial_bound':lower/(NU*t_ref**4),'w_eq':p/rho})
    for cap in (.3,.2976,.4,.107):
        # Scalar-relativistic late-BBN limit, not a likelihood.
        rmax=(cap*7/4)**(.75)/10.75
        result['entropy_targets'].append({'cap':cap,'r_max':rmax,'fraction_retained_max':rmax/r,'xi_initial_max_if_gdark_constant':(rmax/r)**(1/3),'gSM_dec_min_if_gdark_constant':gd/rmax})
    result['high_gsm_106.75_massless_dark_limit']=(4/7)*(8*10.75/106.75)**(4/3)
    s_today=r*S0*HBARC**3
    t,rho,p=min_energy(s_today,86.)
    result['today_bound']={'T_eV':t,'omega_h2_min':rho/HBARC**3/RHO_CRIT,'w_eq':p/rho}
    # If all dark entropy becomes a relativistic equilibrium scalar before relativistic freeze-out:
    # s/n = 2*pi**4/(45*zeta(3)) for one real boson; s_dark today / this gives n.
    zeta3=1.202056903159594
    n0=r*S0/(2*pi**4/(45*zeta3))
    result['relativistic_scalar_only_freezeout']={'n_cm3':n0,'omega_h2':86*n0/RHO_CRIT}
    r_relic=.120*RHO_CRIT/86/S0*(2*pi**4/(45*zeta3))
    result['scalar_relic_entropy_cap']={'r_max':r_relic,'fraction_retained_max':r_relic/r,'xi_initial_max_if_gdark_constant':(r_relic/r)**(1/3)}
    # Additional assumption: total comoving dark count cannot decrease after BBN.
    # Use all 9 degrees as massless bosons to overbound the pressure at z<1.
    s=r*A*10.75
    n=.120*RHO_CRIT/86/S0*A*10.75
    z=exp(4-s/n)
    t=(pi*pi*n/(9*z))**(1/3)
    p=9*t**4/(pi*pi)*z/(1-z)
    result['count_constrained_gate']={'s_over_Tref3':s,'n_max_over_Tref3':n,'s_per_particle_min':s/n,'trial_fugacity':z,'trial_T_over_Tref':t,'rho_lower_over_one_nu':(t*(s+log(z)*n)-p)/NU,'number_destruction':[]}
    for cap in (1.,2.,5.):
        def entropy_upper(n_trial):
            t=cap*NU/(3*n_trial)
            z=pi*pi*n_trial/(9*t**3)
            assert 0<z<1
            p=9*t**4/(pi*pi)*z/(1-z)
            return (cap*NU+p)/t-log(z)*n_trial
        n_min=brentq(lambda nt:entropy_upper(nt)-s,n,.09)
        result['count_constrained_gate']['number_destruction'].append({'rho_cap_in_one_nu':cap,'n_min_over_Tref3':n_min,'n_min_over_n_today_max':n_min/n,'minimum_fraction_destroyed':1-n/n_min})
    return result

if __name__=='__main__':
    print(json.dumps(main(),indent=2))
