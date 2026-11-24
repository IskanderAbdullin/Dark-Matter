"""Exact thermal FD on-shell direct Yukawa decay cold-angle fraction and
independent ultra-relativistic closed-form limit.

Source: χ1 FD(T), isotropic vacuum unpolarised χ1→χ2+real φ.
RESTRICTED F1 only: set Pauli daughter f2=0 and Bose bφ=0, no inverse,
no plasma thermal masses and no FRW/time evolution. This is an initial
direct-production control, NOT a physically frozen DM cold fraction, BEC or Ω.
"""
import math
import json
from scipy.integrate import quad
from scipy.special import expit,spence,zeta
from kinetics.phi_bose_physical_cutcell_weak import (
    M1,M2,MPHI,KSTAR,ESTAR,shell,
)
from kinetics.phi_bose_invariant_parent_support import parent_interval


def ultra_rel_kinematic_fraction(alpha):
    """FD-weighted direct decay fraction kφ/T<α as T/M1→infinity.

    Parent decay count weight q dq/(exp(q)+1), scalar energy is uniform
    over (0,A*p) for isotropic rest-frame decay with
    A=(Eφ*+k*)/M1. Thus direct daughter x=k/T PDF
        (12/(Aπ²))*ln(1+exp(-x/A))
    and CDF=1+(12/π²) Li_2(-exp(-α/A)).
    This distribution is per *direct decay*, not the Bose dark relic.
    """
    if alpha<0:raise ValueError("alpha must be nonnegative")
    A=(ESTAR+KSTAR)/M1
    u=alpha/A
    return float(1.+12/math.pi**2*spence(1.+math.exp(-u)))


def ultra_rel_number_pdf(x):
    """Exact ultra-relativistic conditional daughter energy PDF in x=k/T."""
    if x<0:raise ValueError("x>=0")
    A=(ESTAR+KSTAR)/M1
    return 12/(A*math.pi**2)*math.log1p(math.exp(-x/A))


def ultra_rel_mean_k_over_T():
    """Mean direct scalar x=k/T weighted by parent FD decay rate."""
    A=(ESTAR+KSTAR)/M1
    return float(9*A*zeta(3.,1.)/math.pi**2)


def _fd_decay_weight(q,T):
    """Apart from a positive T,M1,Γ factor, dN/dq of decaying parents."""
    E=math.hypot(q*T,M1)
    return q*q/E*float(expit(-E/T))


def _angle_fraction(q,T,threshold):
    """Exact boosted uniform rest-angle interval, positive massive φ."""
    p=q*T
    if p==0:return float(ESTAR<threshold)
    E,lo,hi,B=shell(p)
    return min(1.,max(0.,(threshold-lo)/(hi-lo)))


def finite_T_direct_cold_fraction(T,alpha):
    """Numerical thermal FD angle+momentum direct-decay fraction (no Bose).

    Adaptive integration starts at exact massive parent p_low(required K).
    Parent FD q tail beyond q_threshold+70 contributes an exponentially
    tiny relative correction. For q_threshold>500, fail rather than
    accidentally report a physically nonzero underflowed source as zero.
    """
    if T<=0 or alpha<=0:raise ValueError("T and alpha positive")
    K=alpha*T
    if K<=MPHI:raise ValueError("Cold energy threshold must exceed mphi")
    bounds=parent_interval(K,MPHI,1e10)
    if bounds is None:raise ValueError("Threshold outside provided support")
    qmin=bounds[0]/T
    if qmin>500:raise ValueError("Boltzmann tail too deep for direct floating integration")
    end=max(75.,qmin+70.)
    N=quad(lambda q:_fd_decay_weight(q,T)*_angle_fraction(q,T,K),
           qmin,end,epsabs=1e-160,epsrel=3e-11,limit=350)[0]
    D=quad(lambda q:_fd_decay_weight(q,T),
           0.,75.,epsabs=1e-12,epsrel=1e-12,limit=350)[0]
    return dict(T_GeV=T,x_threshold=alpha,
                physical_parent_p_min_GeV=bounds[0],
                required_parent_q_min=qmin,
                cold_direct_decay_fraction=N/D,
                parent_q_max_for_numerator=end,
                scope="forward-only direct vacuum isotropic Yukawa, f2=0,bphi=0; not late DM cold fraction")


def finite_T_direct_mean_x(T):
    """Independent physical average scalar energy / mother temperature."""
    if T<=0:raise ValueError("Positive T")
    def source(q):
        return _fd_decay_weight(q,T)
    def heat(q):
        p=q*T
        e,lo,hi,_=shell(p)
        return source(q)*(lo+hi)/(2*T)
    R=quad(source,0.,75.,epsrel=1e-11)[0]
    Q=quad(heat,0.,75.,epsrel=1e-11)[0]
    return Q/R


def audit():
    vals=[]
    for T in (.1,1.,10.,100.,1000.):
        for alpha in (.001,.01,.1):
            vals.append(finite_T_direct_cold_fraction(T,alpha))
    return dict(schema_version=1,
                scope="on-shell FD parent DIRECT unblocked/un-stimulated decay at one fixed T, not full 1↔2 quantum, FRW or relic",
                ultrarelativistic_kinematic_A=(ESTAR+KSTAR)/M1,
                ultrarelativistic_mean_k_over_T=ultra_rel_mean_k_over_T(),
                ultrarelativistic_cold_fractions={str(a):ultra_rel_kinematic_fraction(a)
                                                  for a in (.001,.01,.1)},
                finite_temperature_rows=vals,
                mean_direct_k_over_T={str(T):finite_T_direct_mean_x(T)
                                      for T in (.1,1.,100.,1000.)},
                class_ready=False)


if __name__=="__main__":
    print(json.dumps(audit(),indent=2,allow_nan=False))
