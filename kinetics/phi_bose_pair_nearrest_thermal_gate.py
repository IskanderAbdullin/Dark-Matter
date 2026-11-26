"""Exact near-rest real φ pair threshold and LOG-occupation factor bound.

For equal-mass incoming χ+antiχ ->φ+φ with Dirac m>mφ, to emit
ANY one scalar of lab energy Ephi<=K where mφ<=K<mχ, the incoming
two-particle LAB total energy must obey:
  E_chi+E_anti >=
  2m cosh[arcosh(m/mφ)-arcosh(K/mφ)].
This is the attainable minimum, realized by two perfectly aligned equal
parent momenta with the CM scalar emitted against the boost.

With prescribed FD parents Tchi>0, μ=0 and an equilibrium Bose φ bath
Tphi>0, μphi<mφ, the POSITIVE FORWARD quantum occupation factor of
ANY SUCH EVENT satisfies the pointwise bound
  F=fχ(E1)fanti(E2)(1+bφ(E3))(1+bφ(E4))
   <=exp[-Etotal_min/Tchi] /
      ([1-exp(-(mφ-μφ)/Tφ)]
       [1-exp(-(Etotal_min-K-μφ)/Tφ)]).
We bound only dimensionless quantum occupation factors per event,
NOT integrated invariant collision rate or physical freezeout/BEC.
The reverse factor similarly uses one hard Bose daughter.
At near rest T~1GeV the bound is astronomically tiny; at TeV T
Bose stimulation may be huge, but equilibrium forward/reverse cancel.
"""
import json
import math
from kinetics.phi_bose_physical_cutcell_weak import M1,M2,MPHI

LOG10=math.log(10.)

def minimum_pair_total_energy_for_scalar_below(K,m,mp=MPHI):
    if not (m>mp>0 and mp<=K<m):
        raise ValueError("Require m_parent>energy_K>=m_phi>0")
    zeta=math.acosh(m/mp)
    lower=math.acosh(K/mp)
    return 2*m*math.cosh(max(0.,zeta-lower))

def log10_one_plus_bose(eta):
    if eta<=0:raise ValueError("Subcritical positive Bose η required")
    return -math.log(-math.expm1(-eta))/LOG10

def log10_bose(eta):
    if eta<=0:raise ValueError("Subcritical positive Bose η required")
    # b=exp(-eta)/(1-exp(-eta)).
    return -eta/LOG10+log10_one_plus_bose(eta)

def thermal_pointwise_pair_nearrest_bound(K,m,Tchi,Tphi,mu_phi=0.,
                                          mp=MPHI):
    if not (Tchi>0 and Tphi>0 and mu_phi<mp):
        raise ValueError("Positive thermal T and subcritical real Bose μ")
    Emin=minimum_pair_total_energy_for_scalar_below(K,m,mp)
    high=Emin-K
    if high<=mp:raise ArithmeticError("Other scalar must be on shell")
    eta_low=(mp-mu_phi)/Tphi
    eta_high=(high-mu_phi)/Tphi
    logF=(-Emin/Tchi/LOG10 + log10_one_plus_bose(eta_low)
          +log10_one_plus_bose(eta_high))
    # Reverse scattering occupations R=(1-f1)(1-fanti)b_low b_high,
    # since 0<1-f<=1, and the high Eother>= Emin-K.
    logR=log10_bose(eta_low)+log10_bose(eta_high)
    return dict(parent_mass_GeV=m,scalar_energy_upper_K_GeV=K,
                incoming_pair_total_energy_min_GeV=Emin,
                parent_equal_momentum_min_GeV=math.sqrt((Emin/2)**2-m*m),
                positive_forward_log10_occupation_upper=logF,
                positive_reverse_log10_occupation_upper=logR,
                Tparent_GeV=Tchi,Tphi_GeV=Tphi,mu_phi_GeV=mu_phi,
                scope="exact pointwise thermal FD/BE occupation factor upper bound per on-shell pair event, NOT integrated rate or actual nonthermal Bose burst",
                class_ready=False)

def audit():
    rows=[]
    for m in (M2,M1):
        for K_eV in (86.,86.1,100.,186.):
            for T in (.1,1.,100.,1000.):
                rows.append(thermal_pointwise_pair_nearrest_bound(
                    K_eV*1e-9,m,T,T))
    return dict(schema_version=1,
                scope="pointwise conditional FD/Bose thermal quantum occupation factors and exact near-rest pair kinematic energy floor",
                rows=rows,class_ready=False)

if __name__=="__main__":
    print(json.dumps(audit(),indent=2,allow_nan=False))
