"""Stable 20-MeV chi2 component: conditional symmetric thermal freeze-out gate.

The real-scalar Yukawa y phi(bar chi2 chi1+h.c.) and gauged U(1) vector
q(chi1)=q(chi2) preserve net dark Dirac number U1. In the minimal field
content, the lightest charged fermion chi2 is STABLE even when phi is
the intended 86-eV dark matter. Thus Omega_phi alone is not the total
Omega_dark. This stage critically tests a DEFINED, restricted thermal
branch, NOT the physical five-field freezeout solution.

Tree 2->2 chi2+antichi2->2 identical real phi via virtual chi1:
  sigma*v_Moller = b*v_Moller² + O(v^4,mphi²/mchi²)
  b = y^4 [3 M^6+6 M^5 m+9 M^4 m²+8 M³ m³
           +5 M² m^4+2 M m^5+3 m^6] /
             [96 pi (M²+m²)^4].
This follows by expanding the exact 4x4 Dirac t/u trace at p=0 and
integrating the fully symmetrized, identical-final-phi cross section.
Independent exact spin trace integration tests are external to this module.
The scalar annihilation is P-WAVE (NOT an s-wave constant), and the
extra finite-mass A'->e+e- annihilation has NR s-wave
  a_e=(gD epsilon e)^2 m² /
       [pi ((4m²-mV²)² +mV² GammaV²)].

The conditional toy RD history assumes symmetric chi2/anti chi2 MB
distributions, zero phi chemical potential and phi at THE SAME T,
constant g*rho=g*S=10.75, annihilation-only chi2 rate
<sigma v> = a_e+6 b/x at x=m/T>=3, with no chi1 repopulation,
extra SM channels, nonthermal Bose stimulation, mediator plasma,
neutrino reheating or selfconsistent expanding multi-species physics.
Y=nchi2/s per individual Dirac charge sector:
   dY/dx= -[s/(H x)] <sigma v> (Y²-Yeq²).
Omega_chi2+anti2 h² =2 m Yinf s0/(rho_crit/h²).
Starting from a thermal Yeq at x=3 is an EXPLICIT assumption,
not a demonstrated initial condition. Results are falsification
diagnostics for this branch only, not an allowed/excluded full model.
"""
import math
import json
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
from scipy.interpolate import PchipInterpolator
from functools import lru_cache
from kinetics.phi_bose_pair_independent_vectorized import sigma_independent
from scipy.special import kve

from kinetics.phi_bose_vector_portal_annihilation_resonance import (
    minimal_vector_width,
)
from kinetics.phi_bose_vector_portal_finite_mediator import (
    vector_portal_setup,
)

MCHI2=0.02
MCHI1=0.1
MPHI=86e-9
Y_REF=.001
MPL_NONREDUCED=1.2209e19
S0_CM3=2891.2
RHOC_H2_GEV_CM3=1.05375e-5
H2_DENSITY_CONVERSION=S0_CM3/RHOC_H2_GEV_CM3


def exact_scalar_pwave_b(y=Y_REF,m=MCHI2,M=MCHI1):
    """Exact leading NR scalar chi2+anti2->phi phi coefficient, GeV^-2."""
    if min(m,M)<=0 or y<0:raise ValueError("Positive masses / nonnegative y")
    numerator=(3*M**6+6*M**5*m+9*M**4*m*m+
               8*M**3*m**3+5*M*M*m**4+2*M*m**5+3*m**6)
    return y**4*numerator/(96*math.pi*(M*M+m*m)**4)


def vector_ee_swave_a(gD=.1,epsilon=1e-5,mV=10.,m=MCHI2,GammaV=None):
    """Actual offshell NR s-wave chi2+antichi2->e+e- at s=4m², GeV^-2."""
    if min(m,mV)<=0:return 0.
    if GammaV is None:GammaV=minimal_vector_width(mV,gD,epsilon)["width_GeV"]
    pair=vector_portal_setup(gD,epsilon,mV)["vector_coupling_product"]
    return pair*pair*m*m/(math.pi*((4*m*m-mV*mV)**2+(mV*GammaV)**2))


def mb_yeq_single_dirac(x,gstar_s=10.75,gchi=2):
    """Physical MB particle number / entropy per chi2 charge sector."""
    if min(x,gstar_s,gchi)<=0:raise ValueError("Positive parameters")
    return (45*gchi/(4*math.pi**4*gstar_s)*
            x*x*kve(2,x)*math.exp(-x))


def conditional_symmetric_relic(*,y=Y_REF,m=MCHI2,M=MCHI1,
                                gstar=10.75,gstar_s=10.75,
                                x_initial=3.,x_final=100000.,
                                a_vector=None,rtol=3e-9):
    """Conditional MB single-species Boltzmann ODE, no full SM or five-field."""
    if not (min(m,M,gstar,gstar_s,x_initial)>0 and x_final>x_initial):
        raise ValueError("Positive physical parameters and final x>initial")
    b=exact_scalar_pwave_b(y,m,M)
    a=vector_ee_swave_a(m=m) if a_vector is None else a_vector
    if a<0:raise ValueError("nonnegative vector rate required")
    lam=(2*math.pi**2/45)/1.66*(gstar_s/math.sqrt(gstar))*MPL_NONREDUCED*m
    def rhs(x,Y):
        eq=mb_yeq_single_dirac(x,gstar_s)
        return [-lam/x**2*(a+6*b/x)*(Y[0]**2-eq*eq)]
    sol=solve_ivp(rhs,(x_initial,x_final),
                  [mb_yeq_single_dirac(x_initial,gstar_s)],
                  method="Radau",rtol=rtol,atol=1e-13)
    if not sol.success or sol.y[0,-1]<=0:
        raise ArithmeticError("Nonpositive/unsolved conditional Boltzmann ODE")
    Yfinal=float(sol.y[0,-1])
    return dict(y=y,mchi2_GeV=m,M_virtual_chi1_GeV=M,
                mphi_GeV=MPHI,gstar_rho_assumed=gstar,gstar_S_assumed=gstar_s,
                x_initial=x_initial,x_final=x_final,
                scalar_pwave_b_GeV_minus2=b,
                vector_e_swave_a_GeV_minus2=a,
                is_scalar_s_wave=False,
                Y_chi2_per_charge_final=Yfinal,
                Omega_chi2_plus_anti2_h2=2*m*H2_DENSITY_CONVERSION*Yfinal,
                last_step_rhs_at_final=rhs(x_final,[Yfinal])[0],
                relative_to_OmegaDM_0p120=2*m*H2_DENSITY_CONVERSION*Yfinal/.120,
                physical_thermal_history_solved=False,
                class_ready=False)


def coupling_for_thermal_relic_target(target=.12,low=.01,high=.02):
    """Illustrative y where the *restricted* stable chi2 relic equals target."""
    if not 0<target or not 0<low<high:raise ValueError("positive target bracket")
    func=lambda y:conditional_symmetric_relic(y=y)[
        "Omega_chi2_plus_anti2_h2"]-target
    if func(low)*func(high)>=0:
        raise ValueError("Root not bracketed")
    return float(brentq(func,low,high,xtol=1e-7))



def exact_MB_pair_sigma_v(x,nangle=96):
    """TRUE relativistic MB thermal <sigma vM> from independent vacuum
    t/u Dirac sigma(s) with identical Bose φ final state 1/2 included.

    Shift u=sqrt(s)/T to u=2x+t. The Bessel exponential factors
    e^(-2x) in number-density squared cancel analytically, avoiding
    underflow for x>>1. The integrand still knows exact momenta and
    virtual χ1 pole, rather than substituting 6b/x at x≈7.
    """
    if not (3<=x<=150):
        raise ValueError("Reference intended only for 3<=x<=150")
    T=MCHI2/x
    k2=kve(2,x)
    def f(t):
        u=2*x+t
        return (u*u*(u*u-4*x*x)*
                sigma_independent((T*u)**2,nangle,
                                  y=Y_REF,m_in=MCHI2,m_virtual=MCHI1)*
                kve(1,u)*math.exp(-t))
    integral=quad(f,0,65,epsabs=0,epsrel=3e-7,limit=180)[0]
    return float(integral/(4*x**4*k2*k2))


@lru_cache(maxsize=8)
def _exact_thermal_table(n_grid=32,nangle=96):
    """Fixed logarithmic sample grid, not a runtime-adaptive physical rate."""
    if n_grid<12:raise ValueError("Need >=12 thermal sample points")
    xs=np.geomspace(3.,100.,n_grid)
    sigmas=np.array([exact_MB_pair_sigma_v(float(x),nangle=nangle)
                      for x in xs])
    if not (np.min(sigmas)>0 and np.all(np.isfinite(sigmas))):
        raise ArithmeticError("Nonpositive exact thermal annihilation")
    return xs,sigmas


def conditional_exact_MB_relic(*,y=Y_REF,grid=32,x_initial=3.,
                               gstar=10.75,gstar_s=10.75):
    """Controlled improvement over the leading P-wave truncation at x~7.

    Exact thermal MB χ2χbar2→φφ rate via full vacuum t/u amplitude
    until x=100; from x=100 to infinity use analytic leading P-wave
    tail. The physical next-to-leading tail is O(1/100) relative,
    and the tail itself contributes negligibly to freezeout.
    Other cosmological assumptions match conditional_symmetric_relic.
    """
    if not (x_initial==3. and gstar>0 and gstar_s>0):
        raise ValueError("This exact table is calibrated for x_initial=3")
    xs,sv=_exact_thermal_table(grid)
    interpolation=PchipInterpolator(np.log(xs),np.log(sv))
    lam=((2*math.pi**2/45)/1.66*(gstar_s/math.sqrt(gstar))*
         MPL_NONREDUCED*MCHI2)
    a=vector_ee_swave_a()
    scale=(y/Y_REF)**4
    b=exact_scalar_pwave_b(y)
    def source(x,Y):
        y_eq=mb_yeq_single_dirac(x,gstar_s)
        sigma_v=scale*math.exp(float(interpolation(math.log(x))))+a
        return [-lam/x**2*sigma_v*(Y[0]**2-y_eq*y_eq)]
    sol=solve_ivp(source,(3.,100.),
                  [mb_yeq_single_dirac(3.,gstar_s)],
                  method="Radau",rtol=3e-9,atol=1e-13)
    if not sol.success or sol.y[0,-1]<=0:
        raise ArithmeticError("Failed full t/u thermally averaged freezeout")
    Y100=float(sol.y[0,-1])
    # Conditional no-repopulation integral to infinity. Yeq at x100
    # is exponentially zero; vector a adds the a/100 term.
    Yinf=1/(1/Y100+lam*(a/100.+3*b/(100.**2)))
    omega=2*MCHI2*H2_DENSITY_CONVERSION*Yinf
    return dict(y=y,thermal_grid=grid,x_initial=3.,x_matching_to_NR_tail=100.,
                full_MB_cross_section_reference="independent exact Dirac t/u sigma(s), GG Bessel thermal average",
                exact_thermal_sigma_v_at_x3=scale*float(sv[0]),
                exact_thermal_sigma_v_at_x100=scale*float(sv[-1]),
                scalar_pwave_at_x3=6*b/3.,
                scalar_pwave_at_x100=6*b/100.,
                exact_to_pwave_x3=float(sv[0]/(6*exact_scalar_pwave_b(Y_REF)/3.)),
                exact_to_pwave_x100=float(sv[-1]/(6*exact_scalar_pwave_b(Y_REF)/100.)),
                Y_per_chi2_charge_at_x100=Y100,
                Y_per_chi2_charge_final_with_NR_tail=Yinf,
                Omega_chi2_plus_anti2_h2=omega,
                relative_to_OmegaDM_0p120=omega/.120,
                constant_gstar_toy_only=True,
                full_physical_SM_Bose_FRW_or_relic=False,
                class_ready=False)


def coupling_for_exact_MB_thermal_target(target=.12,low=.01,high=.02,grid=32):
    """Specific artificial y threshold for this restricted thermal branch.

    This is NOT a fit to actual cosmic dark/SM histories or φ production.
    Scaling y changes decay rates and five-field relic, not modeled here.
    """
    if not (target>0 and low>0 and high>low):
        raise ValueError("Positive target and bracket")
    f=lambda y:conditional_exact_MB_relic(y=y,grid=grid)[
        "Omega_chi2_plus_anti2_h2"]-target
    if f(low)*f(high)>=0:raise ValueError("No exact-MB bracket")
    return float(brentq(f,low,high,xtol=1e-7))

def audit():
    rows=[conditional_symmetric_relic(y=y) for y in
          (.001,.002,.004,.008,.01,.012,.014,.02)]
    return dict(schema_version=1,
                goal="Test unavoidable stable Dirac chi2 component in exact minimal U1-conserving real-phi Yukawa plus vectorlike portal",
                all_reaction_rates_NR_truncated=True,
                symmetry="The connected full Yukawa plus vectorlike portal has net dark U1, lightest chi2 cannot decay without additional symmetry-breaking fields",
                target_Omega_DM_h2=.12,
                assumed_gstar_rho_S=10.75,
                rows=rows,
                toy_y_at_Omega_chi2_equal_0p12=coupling_for_thermal_relic_target(),
                caveat="Restricted constant-g* MB symmetric thermal freezeout with p-wave leading scalar and e s-wave offshell, NO physical SM+chi1+phi thermal FRW, no Bose enhancement, no extra mediator channels, no actual halo/LSS likelihood. NOT a full-model exclusion.",
                class_ready=False)


if __name__=="__main__":
    print(json.dumps(audit(),indent=2))
