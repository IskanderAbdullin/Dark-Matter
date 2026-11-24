"""Exact global LAB scalar energy support for equal-mass Dirac pair ->2 real φ.

Incoming χ and antiχ have identical vacuum mass m > mφ, both individual
momentum magnitudes <=P in an isotropic plasma. Any outgoing identical
real scalar energy k satisfies:
  k >= mφ*cosh(max(0, arcosh(m/mφ)-arsinh(P/m))).
This is the ACHIEVABLE global lower boundary, not just a bound.
For P < P*=m*sqrt(m²-mφ²)/mφ, the minimum occurs for aligned equal
parent momenta both P, CM relative momentum zero, backward CM emission:
  k_min=sqrt(P²+m²)-P*sqrt(1-(mφ/m)²).
At P>=P*, one can choose BOTH parents p=P* aligned and produce a
daughter exactly at rest, so k_min=mφ.
Above P*, simply evaluating the endpoint p=P gives a FALSE larger min;
always use the global inverse support formula (max(0,...) rapidities).

Proof: At fixed lab E=sum parent energies and initial invariant s>=4m²,
the minimum daughter lab energy is
  kmin(E,s)=[E-sqrt(E²-s)*sqrt(1-4mφ²/s)]/2.
The term subtracted is maximized at s=2mφ E (if accessible), otherwise
s=4m². For E<=2 sqrt(P²+m²) and P<P*, one has 2mφ E<4m²,
so s=4m² and E=2sqrt(P²+m²) yield the global minimum. The threshold
is achieved by same-velocity aligned parents. When P>=P* a daughter
at rest is accessible with both parents p=P*.

No thermal probability, amplitude/occupation or FRW is calculated.
"""
import json
import math
from kinetics.phi_bose_physical_cutcell_weak import M1,M2,MPHI

def pair_rest_parent_threshold(m,mp=MPHI):
    if not (m>mp>0):raise ValueError("m_parent>m_phi>0 required")
    return m*math.sqrt((m-mp)*(m+mp))/mp

def global_pair_scalar_support_floor(P,m,mp=MPHI):
    if P<0 or not (m>mp>0):
        raise ValueError("P>=0, m_parent>m_phi>0 required")
    rap=math.asinh(P/m)
    daughter=math.acosh(m/mp)
    return mp*math.cosh(max(0.,daughter-rap))

def pair_endpoint_scalar_shell(P,m,mp=MPHI):
    """p1=p2=P aligned, backwards CM. NOT global above Pstar."""
    if P<0 or not (m>mp>0):
        raise ValueError("Physical positive input")
    # Use rapidity difference rather than cancellation at TeV.
    return mp*math.cosh(math.acosh(m/mp)-math.asinh(P/m))

def _generic_min_at_total_lab_energy_and_s(E,s,mp):
    """Independent exact finite two-body lab kinematics at fixed (E,s)."""
    if not (E*E>=s>=4*mp*mp):
        raise ValueError("Physical total E and pair invariant required")
    return (E-math.sqrt(max(0.,E*E-s))*
              math.sqrt(max(0.,1-4*mp*mp/s)))/2

def audit():
    result={}
    for key,m in (("chi2_pair",M2),("chi1_pair",M1)):
        pstar=pair_rest_parent_threshold(m)
        rows=[]
        for P in (0.,2.,8.,16.,100.,1000.,5000.,60000.,130000.):
            floor=global_pair_scalar_support_floor(P,m)
            endpoint=pair_endpoint_scalar_shell(P,m)
            rows.append(dict(P_GeV=P,global_floor_GeV=floor,
                             aligned_endpoint_floor_GeV=endpoint,
                             endpoint_is_global=(P<=pstar)))
        result[key]=dict(parent_mass_GeV=m,pstar_rest_parent_GeV=pstar,
                         rows=rows)
    return dict(schema_version=1,mphi_GeV=MPHI,
                scope="exact physical kinetic global minimum per daughter energy for hard per-incoming-parent P, both equal-mass Dirac and anti; no rate, thermal plasma, BEC or FRW",
                result=result,class_ready=False)

if __name__=="__main__":
    print(json.dumps(audit(),indent=2,allow_nan=False))
