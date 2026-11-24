"""Off-diagonal Dirac Yukawa chi2+anti(chi2)->phi+phi vacuum cross section.

L_int = -y phi bar(chi2) chi1 + h.c., real identical final scalars.
t and u virtual chi1 propagators are combined with PLUS sign, and
initial chi2/anti spin average 1/4, final-state identical scalar 1/2.
Spin trace is evaluated from invariant Minkowski dot products, independent
of any synthetic collision network.

Forward-only Maxwell-Boltzmann thermal reaction-density estimate ignores
inverse phi+phi->chi2+anti(chi2), Bose stimulation, FD corrections,
medium masses, SM/FRW backreaction. NOT a physical relic prediction.
"""
import math
import json
import numpy as np
from scipy.integrate import quad
from scipy.special import kv
from numpy.polynomial.legendre import leggauss

M1=0.1
M2=0.02
MPHI=86e-9
Y=0.001
G_SPIN=2  # chi2 and anti(chi2) separately have two spin states
PLANCK=1.2209e19


def dot(a,b):
    """Mostly-minus Lorentz product of real four-vectors."""
    return float(a[0]*b[0]-np.dot(a[1:],b[1:]))


def scattering_vectors(s,c):
    """Center-of-mass external on-shell four-vectors, real scalar daughters."""
    if s<4*max(M2*M2,MPHI*MPHI):
        raise ValueError("Closed or unphysical pair channel")
    if not (-1<=c<=1):raise ValueError("cos(theta) outside [-1,1]")
    E=math.sqrt(s)/2
    pin=math.sqrt(max(0,E*E-M2*M2))
    pout=math.sqrt(max(0,E*E-MPHI*MPHI))
    transverse=pout*math.sqrt(max(0,1-c*c))
    chi2=np.array([E,0.,0.,pin])
    antichi2=np.array([E,0.,0.,-pin])
    k1=np.array([E,transverse,0.,pout*c])
    k2=np.array([E,-transverse,0.,-pout*c])
    return chi2,antichi2,k1,k2


def trace_quarter(a,q,b,r):
    """1/4 Tr[(slash(a)-m2)(slash(q)+m1)(slash(b)+m2)(slash(r)+m1)].

    From tr(gamma^mu gamma^nu)=4g^munu and trace of four gammas;
    independent direct 4x4 gamma matrix verification is in the tests.
    """
    four=(dot(a,q)*dot(b,r)-dot(a,b)*dot(q,r)+dot(a,r)*dot(q,b))
    two=(M2*M1*dot(a,q)+M1*M1*dot(a,b)+M2*M1*dot(a,r)
         -M2*M1*dot(q,b)-M2*M2*dot(q,r)-M2*M1*dot(b,r))
    return four+two-M2*M2*M1*M1


def spin_averaged_amp2(s,c,y=Y):
    """Unpolarized |M_t+M_u|² for chi2+anti(chi2)->2 real phi."""
    p1,p2,k1,k2=scattering_vectors(s,c)
    q1=p1-k1
    q2=p1-k2
    d1=dot(q1,q1)-M1*M1
    d2=dot(q2,q2)-M1*M1
    qs=(q1,q2)
    ds=(d1,d2)
    total=sum(trace_quarter(p2,q,p1,r)/(dq*dr)
              for q,dq in zip(qs,ds)
              for r,dr in zip(qs,ds))
    out=y**4*total
    if out< -max(1e-29,y**4*1e-13):
        raise ArithmeticError("Spin-averaged squared amplitude negative")
    return max(0.,float(out))


def cross_section(s,y=Y,nangle=36):
    """Total vacuum 2->2 cross section in GeV^-2.

    dσ/dcosθ = |M|²_av/(64πs) *pout/pin,
    INCLUDING 1/2! for two identical real scalar final states.
    """
    if s<=4*max(M2*M2,MPHI*MPHI):return 0.
    E=math.sqrt(s)/2
    pin=math.sqrt(E*E-M2*M2)
    pout=math.sqrt(E*E-MPHI*MPHI)
    x,w=leggauss(nangle)
    amp=np.array([spin_averaged_amp2(s,float(c),y) for c in x])
    return float(np.dot(w,amp)*pout/(64*math.pi*s*pin))


def nchi_MB(T):
    """One chi2 charge sector, spin g=2, ideal unblocked MB particle density."""
    if T<=0:raise ValueError("T>0")
    return G_SPIN*M2*M2*T*kv(2,M2/T)/(2*math.pi**2)


def thermal_forward_rate_MB(T,gstar,y=Y,nangle=36):
    """Forward chi2-anti(chi2)->2phi MB per-chi rate, NOT net collision.

    Gondolo-Gelmini exact MB 2-particle reaction density:
      γ = g² T/(32π⁴) ∫_{4m²}∞ ds σ(s)(s-4m²)√s K1(√s/T).
    σ includes identical final state 1/2, initial spins averaged.
    Rate per χ2 = γ/n_χ2, assuming nχ2=nantiχ2=MB values.
    """
    if T<=0 or gstar<=0:raise ValueError("Positive T and gstar required")
    x=M2/T
    # u=√s/T, so ds=2 T² u du. End u=40 has negligible K1 tail
    # at stated numerical accuracy; narrow high s omitted deliberately.
    def integrand(u):
        s=(T*u)**2
        return (2*T**5*u*u*(u*u-4*x*x)*cross_section(s,y,nangle)
                *kv(1,u))
    moment=quad(integrand,2*x,40.,epsabs=1e-34,
                epsrel=3e-6,limit=140)[0]
    density=(G_SPIN**2*T/(32*math.pi**4))*moment
    n=nchi_MB(T)
    H=1.66*math.sqrt(gstar)*T*T/PLANCK
    return dict(T_GeV=T,gstar=gstar,y=y,MB_particle_density_nchi2_GeV3=n,
                forward_pair_reaction_density_GeV4=density,
                forward_pair_rate_per_chi2_GeV=density/n,
                standard_RD_H_GeV=H,
                forward_pair_rate_over_H=density/n/H,
                class_ready=False)


def audit():
    cs=[dict(sqrt_s_GeV=E,sigma_GeV_minus2=cross_section(E*E))
        for E in (.05,.1,.2,.5,1.)]
    return dict(schema_version=1,
                scope="vacuum spin-averaged off-diagonal real-scalar Yukawa pair amplitude and FORWARD MB thermal rate only, not full quantum net chemical evolution",
                masses_GeV=dict(m1=M1,m2=M2,mphi=MPHI),
                yukawa=Y,
                cross_sections=cs,
                thermal_rows=[thermal_forward_rate_MB(T,g)
                              for T,g in ((.1,17.25),(1.,61.75),(100.,106.75))],
                symmetry_caveat="The pair channel violates isolated decay-only positive number invariants, but conserves the net Dirac U1; exact forward-inverse quantum rate/Γthermal/FRW not computed.",
                class_ready=False)


if __name__=="__main__":
    print(json.dumps(audit(),indent=2,allow_nan=False))
