"""Physical massive Bose first-cell two-moment realizability boundary.

Freeze the *independently normalized physical Yukawa weak source* C at the
thermal initial state and ask when the straight ray M(t)=M(0)+t*C exits the
affine inverse-fugacity DG-P1 family eta(E)=eta_L*hL+eta_R*hR.
A positive non-affine occupation counterexample is then explicitly built.

This is not an evolving self-consistent strong collision solver, condensate,
FRW or cosmological relic prediction.
"""
import math
import json
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from kinetics.phi_bose_exact_cutcell_tangent import tangent
from kinetics.phi_bose_physical_cutcell_weak import MPHI


def _cuts(pmax):
    return sorted({x for x in (0., MPHI, 10*MPHI, 100*MPHI,
                                1.e-6, 1.e-5, 1.e-4, pmax)
                   if 0 <= x <= pmax})


def difference_moments(upper, eta_left, eta_right, temperature=1.):
    """M(eta_left,eta_right)-M(thermal), positive physical p²dp measure.

    Integrate the difference directly, avoiding subtraction of nearly equal
    total thermal moments near a critical endpoint eta_left=0.
    """
    m=MPHI
    assert eta_left >= 0 and eta_right > 0 and upper > m and temperature > 0
    pmax=math.sqrt(upper*upper-m*m)
    cuts=_cuts(pmax)

    def fun(p,i):
        if p==0:
            # For etaL=0: E-m=p²/(E+m), so p²*b -> 2m(U-m)/etaR.
            return (2*m*(upper-m)/eta_right)/(2*math.pi**2) if eta_left==0 and i==0 else 0.
        e=math.hypot(p,m)
        x=p*p/(e+m)  # stable E-m near massive threshold
        h2=x/(upper-m)
        h1=1-h2
        eta=eta_left*h1+eta_right*h2
        delta=1/math.expm1(eta)-1/math.expm1(e/temperature)
        return p*p/(2*math.pi**2)*(h1 if i==0 else h2)*delta

    return np.array([sum(quad(lambda p:fun(p,i),a,b,epsabs=5.e-22,
                              epsrel=8.e-12,limit=300)[0]
                         for a,b in zip(cuts[:-1],cuts[1:])) for i in (0,1)])


def boundary_crossing(n,pmax=2.):
    """The *exact frozen-source* ray intersection with eta_L=0, if approached."""
    t=tangent(n,pmax)
    if not t['first_cell_kinematically_accessible']:
        return dict(n=n,pmax=pmax,support=False,boundary_time=None)
    U=t['first_bin_top_GeV']
    C=np.array(t['source_hat_left_right'])
    slope=np.array(t['deta_left_right'])
    if slope[0]>=0 or C[0]>=0 or C[1]>=0:
        return dict(n=n,pmax=pmax,support=True,
                    Rphi=t['Rphi_first_GeV4'],deta_left=float(slope[0]),
                    boundary_time=None,classification="no left boundary approached locally")

    def diff(er):return difference_moments(U,0.,er)
    def determinant(er):
        d=diff(er)
        return d[0]/C[0]-d[1]/C[1]
    lo=U
    hi=U+max(4*abs(float(slope[1]))*
             t['linear_tangent_positive_eta_limit_GeV_inverse'],2e-6)
    if determinant(lo)*determinant(hi)>=0:
        raise ValueError("Critical face not bracketed; inspect state and test assumptions")
    er=brentq(determinant,lo,hi,xtol=8e-17,rtol=4e-14)
    d=diff(er)
    times=d/C
    tstar=float(np.mean(times))
    return dict(n=n,pmax=pmax,support=True,
                Rphi=t['Rphi_first_GeV4'],
                deta_left=float(slope[0]),deta_right=float(slope[1]),
                upper=U,eta_right_at_left_zero=er,boundary_time=tstar,
                linear_time=t['linear_tangent_positive_eta_limit_GeV_inverse'],
                time_ratio_to_linear=tstar/t['linear_tangent_positive_eta_limit_GeV_inverse'],
                two_moment_relative_crossing_residual=float(abs(times[0]-times[1])/tstar),
                physical_eta_left_at_initial=MPHI,
                C_left_right=C.tolist(),delta_M_left_right=d.tolist())


def dual_gap_at_fixed_right(n,delta_t,pmax=2.):
    """Target ML minus maximum ML in positive *affine-eta* family at fixed MR.

    At fixed MR, dML/detaL = -(G_LL-G_LR²/G_RR)<0 (positive Hessian).
    Therefore etaL=0 maximizes ML, and a positive gap is an exclusion proof.
    """
    a=tangent(n,pmax)
    U=a['first_bin_top_GeV']
    C=np.array(a['source_hat_left_right'])
    target_right=delta_t*C[1]

    def f(er):
        return difference_moments(U,0.,er)[1]-target_right

    lo=U*1e-5
    hi=U+max(1e-5,abs(a['deta_left_right'][1])*delta_t*10)
    if f(lo)*f(hi)>0:
        raise ValueError("Could not bracket eta_R at prescribed target M_R")
    er=brentq(f,lo,hi,xtol=1e-16,rtol=1e-13)
    d=difference_moments(U,0.,er)
    gap=float(delta_t*C[0]-d[0])
    return dict(n=n,dt=delta_t,eta_right_boundary=er,
                target_minus_max_left_M=gap,
                target_right_residual=float(d[1]-target_right),
                classification="OUTSIDE affine Bose family" if gap>0 else "inside-or-boundary")


def positive_constructive_update(n,delta_t,pmax=2.):
    """Explicit admissible non-affine b_new(E)=b_T(E)*(1+aL*hL+aR*hR).

    The physical DOS matrix A_ab=int p²dp/(2pi²)*b_T*h_a*h_b is
    positive definite. Solve A a = delta_t*C to match both moments.
    Endpoint multipliers >0 imply b_new(E)>0 for every E in the cell.
    """
    t=tangent(n,pmax)
    U=t['first_bin_top_GeV']
    m=MPHI
    pmax_bin=math.sqrt(U*U-m*m)
    C=np.array(t['source_hat_left_right'])
    cuts=_cuts(pmax_bin)

    def fun(p,i,j=None):
        if p==0:return 0.
        e=math.hypot(p,m)
        h2=(p*p/(e+m))/(U-m)
        h=(1-h2,h2)
        bose=1/math.expm1(e)
        return p*p/(2*math.pi**2)*bose*h[i]*(h[j] if j is not None else 1)

    A=np.array([[sum(quad(lambda p:fun(p,i,j),l,r,
                         epsabs=2e-20,epsrel=4e-12)[0]
                       for l,r in zip(cuts[:-1],cuts[1:]))
                 for j in (0,1)] for i in (0,1)])
    M=np.array([sum(quad(lambda p:fun(p,i),l,r,
                       epsabs=2e-20,epsrel=4e-12)[0]
                     for l,r in zip(cuts[:-1],cuts[1:]))
                 for i in (0,1)])
    coeff=np.linalg.solve(A,delta_t*C)
    new=M+A@coeff
    return dict(n=n,dt=delta_t,aL_aR=coeff.tolist(),
                multipliers_at_edges=(1+coeff).tolist(),
                initial_M=M.tolist(),target_M=(M+delta_t*C).tolist(),
                constructed_M=new.tolist(),
                max_moment_residual=float(np.max(np.abs(new-M-delta_t*C))),
                positive=bool(np.min(1+coeff)>0))


def audit():
    rows=[]
    for n in (24,32,40):
        c=boundary_crossing(n)
        t=c["boundary_time"]
        rows.append(dict(crossing=c,
                         gap_before=dual_gap_at_fixed_right(n,.8*t),
                         gap_after=dual_gap_at_fixed_right(n,1.1*t),
                         positive_non_affine_at_2x=positive_constructive_update(n,2*t)))
    controls=[boundary_crossing(n,p) for p in (2.,8.,16.)
              for n in (16,24,32,40,48,64)]
    return dict(schema_version=1,model="fixed physical thermal Yukawa WEAK source",
                interpretation="moment-family realizability test only; not self-consistent strong time or BEC",
                rows=rows,cutoff_controls=controls,class_ready=False)


if __name__=="__main__":
    print(json.dumps(audit(),indent=2,allow_nan=False))
