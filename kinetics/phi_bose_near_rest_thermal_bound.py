"""Rigorous thermal 86-eV near-rest Yukawa scalar event rate UPPER bound.

Prescribed FD1/FD2 and strictly subcritical BE phi, NOT a computed rate,
nonthermal Bose burst, BEC, FRW, observational relic, or full strong solver.
"""
import math
import json
from kinetics.phi_bose_physical_cutcell_weak import M1,M2,MPHI,G1,WIDTH
from kinetics.phi_bose_invariant_parent_support import parent_interval


def stable_logsum10(a,b):
    h,l=sorted((a,b),reverse=True)
    return h+math.log10(1+10**(l-h))


def thermal_IR_upper_bound(kmax,temperature_parent=.9,temperature_fermion2=1.,
                            temperature_scalar=1.,mu1=0.,mu2=0.,muphi=0.):
    """Absolute rate |∫_{mφ<=k<=kmax} dW(F−R)| <= 10^log10_R_upper.

    The union of on-shell parent momenta for [mφ,kmax] is the physically
    finite interval from globally two-sided Lorentz inverse rapidity.
    No parent UV cutoff is needed for this bound.
    """
    assert MPHI<kmax<.048
    assert temperature_parent>0 and temperature_fermion2>0 and temperature_scalar>0
    assert muphi<MPHI
    a,b=parent_interval(kmax,MPHI,1e10) # nonbinding numerical cutoff: b<1e7 in these trials
    # Parent FD1 ≤ exp(-(E1-mu1)/T1), E1≥p≥a.
    # Daughter FD2 ≤ exp(-(E2-mu2)/T2), E2=E1-k≥a-kmax.
    logf1=-(a-mu1)/temperature_parent/math.log(10)
    logf2=-(a-kmax-mu2)/temperature_fermion2/math.log(10)
    bmax=1/math.expm1((MPHI-muphi)/temperature_scalar)
    # |J|≤F+R≤f1(1+bmax)+f2*bmax.
    log_flux=stable_logsum10(logf1+math.log10(1+bmax),
                             logf2+math.log10(bmax))
    # ∫dc/2<=1 and ∫p²/E1 dp<=∫p dp=(b²-a²)/2.
    pref=G1*WIDTH*M1*(b*b-a*a)/(4*math.pi**2)
    lograte=math.log10(pref)+log_flux
    return dict(phi_interval_GeV=[MPHI,kmax],
                parent_allowed_union_GeV=[a,b],
                bose_thermal_max=bmax,
                log10_forward_plus_reverse_bound_flux=log_flux,
                prefactor_GeV4=pref,
                log10_absolute_rate_bound_GeV4=lograte,
                log10_absolute_scalar_heat_bound_GeV5=lograte+math.log10(kmax),
                interpretation='rigorous upper bound for prescribed FD1 FD2 BEphi thermal IR only; not nonlinear BEC/FRW')


def audit():
    return dict(schema_version=1,
                rows=[thermal_IR_upper_bound(MPHI+x) for x in
                      (1e-10,1e-9,1e-8,1e-7)],
                class_ready=False)


if __name__=='__main__': print(json.dumps(audit(),indent=2))
