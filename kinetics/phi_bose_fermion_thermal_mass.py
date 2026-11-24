"""One-loop STATIC fermion thermal curvature of the two-Dirac off-diagonal Yukawa.

For L_int=-y phi (bar chi2 chi1+bar chi1 chi2), real phi, and Dirac
mass matrix M(phi)=[[m1,y phi],[y phi,m2]], the one-loop *equilibrium*
fermion free-energy potential is
  V_F(phi,T)=-4T sum_eigen ∫d³p/(2π)^3 log[1+exp(-E_eigen/T)].

At phi=0 its EXACT one-loop static curvature is
  delta m_phi,F^2 = 8y²[m1 I_F(m1,T)-m2 I_F(m2,T)]/(m1-m2),
  I_F=∫d³p/(2π)^3 f_F(E)/E.
Its high-T asymptote for TWO thermal Dirac eigenmodes is y² T²/3.
This STATIC curvature need not equal a thermal quasiparticle POLE mass;
it is not a valid standalone recipe for modifying the on-shell 1↔2 kernel.
All results conditional on two thermalized Dirac species with mu=0 and
vacuum masses. A real scalar quartic lambda*phi^4/4! contributes
lambda*T²/24 in its high-T, one-loop, weakly coupled limit.
"""
import json
import math
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar
from scipy.special import expit

M1=.1
M2=.02
MPHI=86e-9
Y=.001

def fermi_I(m,T):
    """Physical I_F=∫ d³p/(2π)^3 fFD(sqrt(p²+m²))/sqrt(p²+m²)."""
    if m<0 or T<=0:
        raise ValueError("m>=0, T>0")
    x=m/T
    def kernel(q):
        E=math.hypot(q,x)
        return q*q/E*float(expit(-E))
    integral=quad(kernel,0.,46.,epsabs=1e-13,epsrel=2e-12)[0]
    return T*T*integral/(2*math.pi**2)

def fermion_static_curvature(T,y=Y):
    """One-loop static second derivative of thermal fermion free energy at phi=0."""
    if T<=0:raise ValueError("T>0")
    return 8*y*y/(M1-M2)*(M1*fermi_I(M1,T)-M2*fermi_I(M2,T))

def fermion_mass_eigenvalues(phi,y=Y):
    """Real eigenmasses at scalar background (only small phi used in audit)."""
    a=np.linalg.eigvalsh(np.array([[M1,y*phi],[y*phi,M2]],float))
    return a

def fermion_thermal_potential(phi,T,y=Y):
    """Independent direct two-Dirac fermion determinant at real scalar phi."""
    if T<=0:raise ValueError("T>0")
    eig=fermion_mass_eigenvalues(phi,y)
    def per_mode(m):
        x=abs(m)/T
        return quad(lambda q:q*q*math.log1p(
                       math.exp(-math.hypot(q,x))),
                    0.,46.,epsabs=2e-13,epsrel=3e-12)[0]
    return -4*T**4/(2*math.pi**2)*sum(per_mode(m) for m in eig)

def total_static_curvature(T,lam=0.,y=Y):
    """Approximate static curvature with high-T scalar quartic tadpole."""
    if lam<0:raise ValueError("lam must be nonnegative for this bounded scalar model")
    return MPHI**2+fermion_static_curvature(T,y)+lam*T*T/24.

def quartic_needed_for_nonnegative_static_curvature(T,y=Y):
    """Minimum lambda in scalar V=lambda*phi^4/4! at this T, under assumptions."""
    return max(0.,-24*(MPHI**2+fermion_static_curvature(T,y))/(T*T))

def critical_curvature_roots(y=Y,lam=0.):
    """Illustrative two crossings of static curvature for benchmark masses.

    Brackets are deliberately restricted to the tested zero-mu thermal
    Dirac regime, not cosmological phase transition temperatures.
    """
    return (
        float(brentq(lambda t:total_static_curvature(t,lam,y),
                     .002,.005,xtol=1e-14)),
        float(brentq(lambda t:total_static_curvature(t,lam,y),
                     .03,.04,xtol=1e-14)),
    )

def critical_quartic_envelope(y=Y):
    """Worst required lambda over 0.004<=T<=0.035GeV; scanned interval only."""
    z=minimize_scalar(
        lambda T: -quartic_needed_for_nonnegative_static_curvature(T,y),
        bounds=(.004,.035),method="bounded",options={"xatol":1e-12})
    return dict(temperature_GeV=float(z.x),
                min_quartic_lambda_phi4_over_4factorial=float(-z.fun),
                bounds_GeV=[.004,.035],scalar_tadpole="lambda*T²/24")

def naive_phi_only_highT_decay_closure(y=Y):
    """Only hypothetical phi-only m_eff²=m_vac²+y² T²/3 equality to m1-m2.

    Does NOT account for fermion thermal mass, Landau cuts, 2to2, plasma
    dispersion, finite SM heat or actual in-medium decay thresholds.
    """
    return math.sqrt(3*((M1-M2)**2-MPHI**2))/y

def audit():
    points=[]
    for T in (.002,.003,.005,.01,.018598494875,.02,.03,
              .035,.04,.1,1.,10.,100.):
        d=fermion_static_curvature(T)
        tot=total_static_curvature(T)
        points.append(dict(T_GeV=T,fermion_static_delta_m2_GeV2=d,
                           highT_y2T2_over3_GeV2=Y**2*T*T/3,
                           ratio_to_highT=d/(Y**2*T*T/3),
                           total_static_curvature_GeV2=tot,
                           signed_curvature="negative" if tot<0 else "positive",
                           static_mass_if_positive_GeV=math.sqrt(tot) if tot>0 else None,
                           required_lambda_phi4_over_4factorial=
                               quartic_needed_for_nonnegative_static_curvature(T)))
    return dict(schema_version=1,
                model="two Dirac mass-eigenstates M(phi)=[[m1,yphi],[yphi,m2]], real scalar phi, common thermal FD mu=0, vacuum m1=.1 m2=.02 GeV, mphi=86eV y=.001",
                formal_status="one-loop STATIC effective potential conditional on equilibrated Dirac bath; not full pole dispersion or physical FRW phase transition",
                highT_asymptote_delta_m2_GeV2="y² T²/3",
                mu_zero_fermion_curvature_roots_T_GeV=
                    critical_curvature_roots(),
                minimum_scalar_quartic_envelope=
                    critical_quartic_envelope(),
                hypothetical_phi_only_decay_closure_T_GeV=
                    naive_phi_only_highT_decay_closure(),
                samples=points,
                references=[
                    "https://doi.org/10.1103/PhysRevD.9.3320",
                    "https://arxiv.org/abs/hep-ph/9901312",
                ],class_ready=False)

if __name__=="__main__":
    print(json.dumps(audit(),indent=2,allow_nan=False))
