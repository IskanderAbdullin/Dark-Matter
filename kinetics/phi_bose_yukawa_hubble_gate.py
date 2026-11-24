"""Physical Yukawa tree width versus radiation-era Hubble time: consistency gate.

For L_int=-y φ \barχ2 χ1+h.c. with Dirac fermion χ1→χ2+real φ,
unpolarized tree average |M|²=y²[(M1+M2)²−mφ²], and
Γrest=|M|²sqrt(λ)/(16π M1³).

FRW H=1.66 sqrt(g*ρ) T²/M_Pl is a standard RD approximation at fixed g*ρ.
The MB average time-dilated decay rate Γrest K1(M1/T)/K2(M1/T)
is an independent *one-way test-particle* time-scale comparison; full quantum
inverse reactions, thermal masses, SM portal and g*(T) are NOT solved.
"""
import math
import json
from scipy.optimize import brentq
from scipy.special import kve

M1=0.1
M2=0.02
MPHI=86e-9
YBENCH=1e-3
MPLANCK_GEV=1.2209e19
HBAR_GEV_S=6.582119569e-25


def rest_width(y=YBENCH):
    """Yukawa vacuum rest width for spin-averaged Dirac→Dirac+real scalar."""
    if y<0:raise ValueError("y must be nonnegative")
    if M1<=M2+MPHI:raise ValueError("Decay not kinematically allowed")
    lamb=(M1*M1-(M2+MPHI)**2)*(M1*M1-(M2-MPHI)**2)
    spin_average_m2=y*y*((M1+M2)**2-MPHI*MPHI)
    return spin_average_m2*math.sqrt(lamb)/(16*math.pi*M1**3)


def hubble_rd(T,gstar):
    """Natural-unit H (GeV) in radiation domination with nonreduced M_Pl."""
    if T<=0 or gstar<=0:raise ValueError("T,gstar positive")
    return 1.66*math.sqrt(gstar)*T*T/MPLANCK_GEV


def mb_mean_rest_time_dilation(T):
    """Exact ideal MB thermal M1/E1 average K1(x)/K2(x)."""
    if T<=0:raise ValueError("Positive T required")
    x=M1/T
    return float(kve(1,x)/kve(2,x))


def rate_gate(T,gstar,y=YBENCH):
    """Diagnostic only, not a full thermal reaction-equilibration condition."""
    gamma=rest_width(y)
    H=hubble_rd(T,gstar)
    dilation=mb_mean_rest_time_dilation(T)
    return dict(T_GeV=T,gstar_rho=gstar,y=y,
                Gamma_rest_GeV=gamma,H_GeV=H,
                Gamma_rest_over_H=gamma/H,
                MB_averaged_rest_time_dilation=dilation,
                Gamma_MB_per_parent_over_H=gamma*dilation/H,
                benchmark_y_for_kappa_rest_0p1=
                    y*math.sqrt(.1*H/gamma) if y>0 else None,
                benchmark_y_for_MB_mean_ratio_0p1=
                    y*math.sqrt(.1*H/(gamma*dilation)) if y>0 else None)


def mb_crossover_T(gstar=106.75,y=YBENCH):
    """Temperature where thermal test-parent MB decay rate equals RD H.

    gstar is frozen to supplied constant; no claim of full equilibrium
    or actual cosmological decoupling temperature.
    """
    def f(T):
        return rest_width(y)*mb_mean_rest_time_dilation(T)/hubble_rd(T,gstar)-1.
    return float(brentq(f,1.,1e5,rtol=1e-12))


def audit():
    cases=[rate_gate(T,g) for T,g in
           [(0.1,17.25),(0.1,61.75),(1.,61.75),
            (100.,106.75),(500.,106.75),(1000.,106.75)]]
    return dict(schema_version=1,
                scope="physical vacuum spin-averaged tree width and conditional MB test-parent decay/Hubble diagnostic, NOT Γthermal full quantum, FRW Boltzmann relic or validated SM portal",
                masses_GeV=dict(m1=M1,m2=M2,mphi=MPHI),
                benchmark_y=YBENCH,
                benchmark_tau_rest_seconds=HBAR_GEV_S/rest_width(),
                MB_ratio_one_T_GeV_at_fixed_gstar106p75=mb_crossover_T(),
                cases=cases,
                existing_FRW_pilot_kappa_range=[1e-4,.2],
                class_ready=False)

if __name__=="__main__":
    print(json.dumps(audit(),indent=2,allow_nan=False))
