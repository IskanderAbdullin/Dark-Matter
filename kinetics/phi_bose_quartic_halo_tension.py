"""Conditional 86eV real-scalar thermal quartic versus Born cluster scattering.

For V=lambda*phi**4/4!, one real identical-scalar elastic final state,
sigma_NR=lambda**2/(128*pi*mphi**2).
lambda stabilizes STATIC fermion-thermal curvature only under the
previously specified equilibrium two-Dirac mass-mixing assumptions.
Cluster sigma/m reference values are illustrative, NOT likelihood bounds.
"""
import json
import math
from kinetics.phi_bose_fermion_thermal_mass import (
    MPHI,critical_quartic_envelope,total_static_curvature
)
HBARC_GEV_CM=0.1973269804e-13
GEV_MASS_G=1.78266192e-24
CM2_G_PER_GEV_MINUS3=HBARC_GEV_CM**2/GEV_MASS_G

def sigma_contact_GeVminus2(lam,m=MPHI):
    if lam<0 or m<=0:raise ValueError("lambda nonnegative, mass positive")
    return lam**2/(128*math.pi*m*m)

def sigma_over_m_cm2_g(lam,m=MPHI):
    return sigma_contact_GeVminus2(lam,m)*CM2_G_PER_GEV_MINUS3/m

def lambda_for_reference_sigma_m(bound,m=MPHI):
    if bound<=0 or m<=0:raise ValueError("positive inputs required")
    return math.sqrt(bound*128*math.pi*m**3/CM2_G_PER_GEV_MINUS3)

def fixed_lambda_mass_for_reference_sigma_m(lam,bound):
    if lam<0 or bound<=0:raise ValueError("positive inputs required")
    return (lam**2*CM2_G_PER_GEV_MINUS3/(128*math.pi*bound))**(1/3)

def audit():
    gate=critical_quartic_envelope()
    lam=gate["min_quartic_lambda_phi4_over_4factorial"]
    refs=[]
    for limit in (.2,1.,2.,5.):
        upper=lambda_for_reference_sigma_m(limit)
        refs.append(dict(illustrative_bound_cm2_g=limit,
                         lambda_max_under_reference=upper,
                         required_over_reference_lambda=lam/upper,
                         curvature_at_20MeV_with_reference_lambda_GeV2=
                            total_static_curvature(.02,lam=upper),
                         fixed_lambda_mass_scaling_eV=
                            fixed_lambda_mass_for_reference_sigma_m(lam,limit)*1e9))
    return dict(schema_version=1,mphi_vacuum_eV=MPHI*1e9,
                scope="conditional equilibrated two-Dirac static one-loop Z2 curvature vs dominant nonrel dilute real-scalar Born contact halo, NOT a full cluster likelihood",
                min_thermal_curvature_lambda=lam,
                sigma_NR_GeVminus2=sigma_contact_GeVminus2(lam),
                sigma_over_mass_cm2_g=sigma_over_m_cm2_g(lam),
                viscosity_sigma_over_mass_cm2_g=2*sigma_over_m_cm2_g(lam)/3,
                reference_scenarios=refs,
                class_ready=False)

if __name__=="__main__":
    print(json.dumps(audit(),indent=2,allow_nan=False))
