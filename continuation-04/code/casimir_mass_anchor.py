"""A FORMAL Casimir high-occupation coefficient match to a homogeneous relic patch.
NOT a measured halo-mass relation or successful full collisionless mapping.
Needs casimir_ir_tail.py output at ../results/casimir_ir_tail.json.
"""
import json
from pathlib import Path
import numpy as np
from scipy.special import zeta
from scipy.constants import G,c,hbar,elementary_charge,parsec
M_sun=1.98847e30
nuc=zeta(1.5)
R=600*parsec
V0=44100
m_rounded=86*elementary_charge/c**2
Tb0_eV=9.7519e-5
p0=Tb0_eV*elementary_charge/c
V_L=3*np.sqrt(3)*np.pi**2*R**3*(m_rounded*V0/p0)**3
R_L=(3*V_L/(4*np.pi))**(1/3)
nb0=zeta(3)/np.pi**2*(p0/hbar)**3
M_L=nb0*m_rounded*V_L
M0=V0**2*R/G
mcrit=(np.sqrt(np.pi/2)*hbar**3/(G*R**2*V0*nuc))**.25
V_L_crit=3*np.sqrt(3)*np.pi**2*R**3*(mcrit*V0/p0)**3
M_L_crit=nb0*mcrit*V_L_crit
ratio=3*np.sqrt(3)*zeta(3)*np.sqrt(np.pi/2)/nuc
B=128/(35*np.sqrt(np.pi)*nuc)

def run():
    result_path=Path(__file__).resolve().parent.parent/"results"/"casimir_ir_tail.json"
    results=json.loads(result_path.read_text())
    rows=results.get("rows") or results["numerical_check"]["values"]
    errors=[]
    for rec in rows:
        if rec["F"]<100:continue
        estimate=1+B*np.sqrt(rec["L"])
        errors.append(dict(F=rec["F"],direct_ratio=rec["ratio"],
           analytic_first_order=estimate,abs_error=rec["ratio"]-estimate))
    assert abs(M_L_crit/M0-ratio)<1e-12
    assert abs(errors[-1]["abs_error"])<1e-3
    return dict(R_kpc=.6,V0_km_s=44.1,
       source_rounded_particle_mass_eV=86,
       exact_critical_particle_mass_eV=mcrit*c*c/elementary_charge,
       primordial_momentum_T_com_eV=Tb0_eV,
       FORMAL_Lagrangian_radius_kpc=R_L/(1000*parsec),
       FORMAL_Lagrangian_mass_Msun=M_L/M_sun,
       FORMAL_exact_critical_patch_mass_Msun=M_L_crit/M_sun,
       mass_scale_V0sq_R_over_G_Msun=M0/M_sun,
       FORMAL_critical_mass_ratio=ratio,
       high_F_next_order_coefficient_B=B,
       correction_check=errors,tests_passed=2,
       warning="Very restrictive fine-grained closed-patch matching; full level-set spectra differ, so NOT a real halo-mass prediction.")

if __name__=="__main__":print(json.dumps(run(),indent=2))
