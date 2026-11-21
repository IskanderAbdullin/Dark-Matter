"""WKB self-consistency of the locally homogeneous Bose–Jeans calculation.
Run: python wkb_consistency.py; requires numpy/mpmath.
This DOES NOT solve the global spectral stability problem.
"""
import json
import numpy as np
from mpmath import mp
mp.dps=40
NUC=mp.zeta(mp.mpf('1.5'))
R_PC=600.
def A(delta):
    z=mp.exp(-delta)
    return mp.polylog(mp.mpf('.5'),z)/mp.polylog(mp.mpf('1.5'),z)
def t_from_nu(nu):
    lo=mp.mpf(0);hi=mp.mpf(100)
    for _ in range(140):
        md=(lo+hi)/2
        if mp.polylog(mp.mpf('1.5'),mp.exp(-md))>mp.mpf(str(nu)):lo=md
        else:hi=md
    return (lo+hi)/2
def half_response_radius(nu):
    d0=t_from_nu(nu)
    A0=A(d0)
    lo=d0;hi=mp.mpf(100)
    for _ in range(140):
        mid=(lo+hi)/2
        if A(mid)>A0/2:lo=mid
        else:hi=mid
    y=(lo+hi)/2-d0
    rh=R_PC*mp.sqrt(6*y)
    lam=2*mp.pi*R_PC/mp.sqrt(A0)
    return dict(nu=nu,delta=float(d0),central_A=float(A0),
        half_response_radius_pc_small_r_potential=float(rh),
        homogeneous_central_Jeans_lambda_pc=float(lam),
        lambda_over_half_response_radius=float(lam/rh))
def main():
    rows=[half_response_radius(nu) for nu in (2.2,2.5,2.6,2.61,2.612)]
    x=rows[3]
    assert abs(x['half_response_radius_pc_small_r_potential']-1.707013)<1e-5
    assert abs(x['lambda_over_half_response_radius']-69.4014)<.001
    d=x['delta'];R=R_PC
    b=2.356077011282519e-13
    xicross=2*np.sqrt(6)*np.sqrt(np.pi)/float(NUC)*b
    rcontact=R*xicross
    lamcontact=2*np.pi*R*np.sqrt(b)
    assert lamcontact/rcontact>1e6
    result=dict(model="WKB applicability of local homogeneous Jeans diagnostic in stratified source Bose halo",
      source_R_pc=R_PC,nu_critical=float(NUC),rows=rows,
      source_nu2_61={**x,
        "alternative_inner_delta_radius_pc":float(R*np.sqrt(6*d)),
        "note":"Susceptibility varies over ~1-2pc, homogeneous marginal wavelength ~118pc."},
      exactly_critical_finite_HF={'b':b,'contact_layer_pc_order':float(rcontact),
       'homogeneous_HF_Jeans_pc':float(lamcontact),
       'ratio_wavelength_to_layer':float(lamcontact/rcontact)},
      tests_passed=3,
      caveat="Central potential y≈xi²/6; not a global spherical stability conclusion.")
    print(json.dumps(result,indent=2))
if __name__=="__main__":main()
