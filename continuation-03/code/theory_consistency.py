"""Independent Bose halo theory checks; NO network, no GitHub Actions.
Requires Python 3.10+, scipy, mpmath.
Use: python theory_consistency.py > theory_check_output.json
Comparisons to the JCAP original TeX formulae; no cosmological/observational inference.
"""
import json, math
from scipy.constants import G, hbar, c, elementary_charge, parsec
from scipy.integrate import quad
from scipy.special import zeta
from mpmath import mp
mp.dps=30
AU=149597870700.0
M_EV=elementary_charge/c**2
def bose_moments(nu):
    zc=mp.zeta(mp.mpf('1.5'))
    if nu>=float(zc):
        z=mp.mpf(1)
    else:
        left,right=mp.mpf(0),mp.mpf(1)
        for i in range(100):
            mid=(left+right)/2
            if mp.polylog(mp.mpf('1.5'),mid)>nu:right=mid
            else:left=mid
        z=(left+right)/2
    l3=mp.polylog(mp.mpf('1.5'),z)
    l5=mp.polylog(mp.mpf('2.5'),z)
    l7=mp.polylog(mp.mpf('3.5'),z)
    return dict(nu=float(nu),z=float(z),
                sigma_over_V0=float(mp.sqrt(l5/l3)),
                excess_kurtosis=float(3*l7*l3/l5**2-3))
def occupied_fraction(F):
    qmax=math.sqrt(math.log1p(1/F))
    def integrand(q):
        return 1-q*q/2+q**4/12 if q<1e-5 else q*q/math.expm1(q*q)
    return quad(integrand,0,qmax,epsabs=1e-13)[0]/(math.sqrt(math.pi)/4*zeta(1.5))
def mixed_growth(fw,z_eq=3400,z_end=5):
    p=(math.sqrt(25-24*fw)-1)/4
    S=(1-fw)*((1+z_eq)/(1+z_end))**(p-1)
    return dict(fw=fw,growth_exponent=p,conditional_power_ratio=S*S)
def main():
    m=86*M_EV;R=600*parsec;V0=44100.;sigma_over_m=.125
    rho=V0*V0/(4*math.pi*G*R*R);n=rho/m
    lam=math.sqrt(2*math.pi)*hbar/(m*V0)
    scattering_length=math.sqrt((sigma_over_m*m)/(8*math.pi))
    radius_firstzero=math.pi*math.sqrt(hbar*hbar*scattering_length/(G*m**3))
    radius_lateformula=(sigma_over_m*hbar**4/(8*math.pi*G**2*m**5))**.25
    g=4*math.pi*hbar*hbar*scattering_length/m
    b=2*g*n/(m*V0*V0)
    jeans_fluid=2*math.pi*R*math.sqrt(b)
    nuc=zeta(1.5);nu=n*lam**3
    result=dict(
      benchmark="m=86eV, F583-4 R=0.6kpc V0=44.1km/s, sigma/m=1.25 cm2/g chosen illustrative",
      first_zero_radius_AU=radius_firstzero/AU,
      source_late_radius_AU=radius_lateformula/AU,
      radius_ratio=radius_firstzero/radius_lateformula,
      thermal_debroglie_m=lam,interparticle_separation_m=n**(-1/3),
      phase_space_n_lambda3=nu,zeta3over2=nuc,
      mass_min_at_fixed_R_V0_eV=86*(nu/nuc)**.25,
      critical_velocity_moments=bose_moments(nuc),
      classical_low_nu_moments=bose_moments(.03),
      critical_occupancy_fractions={str(F):occupied_fraction(F) for F in [1,3,10,30,100,1000]},
      hartree_fock_b=b,
      isothermal_sound_speed_ms=V0*math.sqrt(b),
      isothermal_fluid_Jeans_pc=jeans_fluid/parsec,
      Jeans_to_TF_first_zero=jeans_fluid/radius_firstzero,
      conditional_smooth_warm_mixed_growth=[mixed_growth(x) for x in [.01,.03,.05,.1,.2,.5]],
      limitation="illustrative mathematical checks only, not 86eV exclusion or galaxy fit"
    )
    assert abs(result['radius_ratio']-math.pi)<1e-12
    assert abs(result['Jeans_to_TF_first_zero']-math.sqrt(8))<1e-10
    assert abs(result['critical_velocity_moments']['sigma_over_V0']-.716597827791)<1e-10
    assert abs(result['critical_velocity_moments']['excess_kurtosis']-1.906875927938)<1e-10
    assert abs(result['critical_occupancy_fractions']['10']-.2625014844)<1e-8
    assert abs(result['mass_min_at_fixed_R_V0_eV']-86.2049585235)<1e-6
    print(json.dumps(result,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
