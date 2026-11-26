"""Collisionless thermal PSD diagnostic, not a perturbation likelihood."""
import json
from math import pi, sqrt, expm1, exp, asinh
from scipy.integrate import quad
from scipy.special import zeta, kve
HBARC=1.973269804e-5 # eV cm
RHOC=1.05375e4 # eV cm^-3
OMEGA_H2=.120
H0=67.4 # km/s/Mpc
C=299792.458 # km/s
OM=.315
H=H0/100
OR=2.4728e-5*(1+.22710731766*3.044)/H**2
OL=1-OM-OR
AEQ=OR/OM


def moment(k,fd):
    return float((1-2**(-k) if fd else 1)*zeta(k+1)*__import__('math').factorial(k))


def weight(x,fd):
    return x*x/(exp(x)+1 if fd else expm1(x))


def species(m,g,fd,fraction=1.):
    n=OMEGA_H2*RHOC*fraction/m
    t=(2*pi*pi*n*HBARC**3/(g*moment(2,fd)))**(1/3)
    r=t/m
    norm=moment(2,fd)
    def average(func):
        return quad(lambda x:weight(x,fd)*func(x),0,60,epsabs=1e-12,epsrel=1e-9)[0]/norm
    # Comoving horizon for momentum q, assuming standard radiation + matter + Lambda background.
    # dchi/da = v/(a^2 H(a)); a^2 H = H0 sqrt(Or+Om*a+Ol*a^4).
    def horizon(a):
        def f(x):
            u=r*x
            return C/H0*quad(lambda aa:u/sqrt(u*u+aa*aa)/sqrt(OR+OM*aa+OL*aa**4),0,a,epsabs=1e-9,epsrel=1e-8,points=[min(u,a)])[0]
        return average(f)
    def vmean(a):return C*average(lambda x:r*x/sqrt((r*x)**2+a*a))
    def vrms(a):return C*sqrt(average(lambda x:(r*x)**2/((r*x)**2+a*a)))
    rd_horizon= C/H0/sqrt(OR)*average(lambda x:r*x*asinh(AEQ/(r*x)))
    return {'mass_eV':m,'g':g,'fd':fd,'fraction':fraction,'n_cm3':n,'T0_eV':t,'q_rms_eV':t*sqrt(moment(4,fd)/norm),'v_rms_today_km_s':vrms(1.),'v_mean_eq_km_s':vmean(AEQ),'v_rms_eq_km_s':vrms(AEQ),'free_stream_to_eq_Mpc':horizon(AEQ),'free_stream_to_z3_Mpc':horizon(.25),'pure_RD_to_eq_Mpc':rd_horizon}


def main():
    be=species(86.,1,False)
    fd=species(5700.,2,True)
    # Chosen diagnostic, not observationally inferred: match NR q/m rms to 5.7keV FD all-DM case.
    q_target=86*fd['q_rms_eV']/5700
    n=be['n_cm3']*HBARC**3
    # Markov: at least half the particles in q<=sqrt(2)*q_rms. Occupation average over that physical phase-space ball.
    q_ball=sqrt(2)*q_target
    f_avg_min=.5*n*6*pi*pi/q_ball**3
    ratio=be['q_rms_eV']/86/(fd['q_rms_eV']/5700)
    result={'H0_km_s_Mpc':H0,'Omega_r':OR,'Omega_m':OM,'a_eq':AEQ,'scalar_BE':be,'reference_FD':fd,'rms_q_over_mass_ratio':ratio,'diagnostic_q_target_eV':q_target,'occupation_ball_average_lower_if_matching_rms':f_avg_min,'cold_plus_BE_rms_match_fraction_diagnostic':ratio**(-6/5),'mixed_tail_examples':[species(86.,1,False,f) for f in (.1,.01,.001)]}
    # Variational entropy upper bound at fixed second moment, arbitrary Bose PSD.
    # log Z/V=A*beta^(-3/2) for trial energy q^2 and zero chemical potential.
    m2=n*q_target*q_target
    aa=float(zeta(2.5)/(8*pi**1.5))
    beta=(1.5*aa/m2)**.4
    smax=2.5*aa*beta**(-1.5)
    nex=float(zeta(1.5)/(8*pi**1.5)*beta**(-1.5))
    # Upper bound on entropy in heavy particles at <=1% of target energy.
    # Fixed trial T=20 MeV/90, using physical FD pressure log integral.
    t_h=20e6/90
    p_h=0.
    from math import log1p
    for mass in (20e6,100e6):
        z=mass/t_h
        p_h+=4*t_h**4/(2*pi*pi)*quad(lambda x:x*x*log1p(exp(-sqrt(x*x+z*z))),0,100,epsabs=1e-100,epsrel=1e-9)[0]
    p_h_upper=sum(4*mass*mass*t_h*t_h*kve(2,mass/t_h)*exp(-mass/t_h)/(2*pi*pi) for mass in (20e6,100e6))
    assert p_h<=p_h_upper*(1+1e-8)
    s_h=(.01*OMEGA_H2*RHOC*HBARC**3+p_h_upper)/t_h
    s_init=(7.991670795139369/61.75)*2891.2
    result['cold_entropy_gate']={'trial_beta_eV_minus2':beta,'q2_density_eV5':m2,'scalar_entropy_upper_cm3':smax/HBARC**3,'saturating_excited_count_cm3':nex/HBARC**3,'saturating_cold_count_fraction':1-nex/n,'heavy_entropy_upper_cm3_for_1pct_energy':s_h/HBARC**3,'total_entropy_upper_cm3':(smax+s_h)/HBARC**3,'retained_fraction_max_from_initial_thermal_benchmark':(smax+s_h)/HBARC**3/s_init}
    return result

if __name__=='__main__':print(json.dumps(main(),indent=2))
