"""Repulsive-TF Bose star size, single-orbit host tidal floor and extended-lens scale.
These are EXAMPLES, not a predicted function of masses or lensing event rates.
"""
import json
import numpy as np
from scipy.constants import G,hbar,c,elementary_charge,parsec
M_sun=1.98847e30
AU=149597870700.
year=365.25*86400
m=86*elementary_charge/c**2
sigma_per_m=.125
scattering_a=np.sqrt(sigma_per_m*m/(8*np.pi))
R=np.pi*np.sqrt(hbar*hbar*scattering_a/(G*m**3))
M_TF_scale=m*R/scattering_a

def M_tidal(Vc_km_s=44.1,orbit_kpc=1):
    omega=Vc_km_s*1e3/(orbit_kpc*1e3*parsec)
    return 2*omega**2*R**3/G

def Einstein_R(M,D_eff_kpc=8):
    return np.sqrt(4*G*M*D_eff_kpc*1e3*parsec/c**2)

def t_dyn_yr(M):return np.sqrt(R**3/(G*M))/year

def run():
    Mt=M_tidal()
    Mdyn10=R**3/(G*(10e9*year)**2)
    M_point=c*c*R**2/(4*G*8e3*parsec)
    rows=[]
    for M in (1e15,1e18,1e20,1e21,1e23,M_sun,100*M_sun,300*M_sun):
        rows.append(dict(M_kg=float(M),M_Msun=float(M/M_sun),
           t_dyn_yr=float(t_dyn_yr(M)),
           R_Einstein_AU=float(Einstein_R(M)/AU),
           R_star_over_R_Einstein=float(R/Einstein_R(M))))
    assert abs(R/AU-133.4461348418)<1e-7
    assert Mt>Mdyn10
    assert R/Einstein_R(M_sun)>10
    return dict(m_eV=86,sigma_over_m_cm2_g=1.25,
        radius_first_zero_AU=R/AU,
        scattering_length_m=scattering_a,
        FORMAL_TF_Na_over_R_eq1_mass_kg=M_TF_scale,
        chosen_host=dict(Vc_km_s=44.1,r_kpc=1,
                         orbit="circular, locally flat VC, spherical host"),
        conditional_tidal_min_kg=Mt,
        conditional_tidal_min_Msun=Mt/M_sun,
        dynamical_10Gyr_min_kg=Mdyn10,
        ratio_tidal_to_formal_TF_scale=Mt/M_TF_scale,
        lensing_D_eff_kpc=8,lensing_point_mass_crossover_Msun=M_point/M_sun,
        trial_masses=rows,tests_passed=3,
        caveat="Not general for eccentric orbits, true host potential, changing star radius, actual Bose-star formation or mass function.")
if __name__=="__main__":print(json.dumps(run(),indent=2))
