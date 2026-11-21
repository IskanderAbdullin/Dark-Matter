"""Exactly critical lowered Bose family across eta=psi0.
Requires scipy/numpy and sibling critical_lowered.py.
No observational fitting or dynamical stability proof.
"""
import json
import numpy as np
from scipy.special import zeta
from critical_lowered import solve
def run():
    rows=[]
    for eta in (.5,1,2,3,4,5,6):
        Ic,rt,Qt,_,_=solve(eta,nq=240)
        rows.append(dict(eta=eta,Icritical=Ic,
          linear_center_slope=np.pi/(2*np.sqrt(6)*Ic),
          central_phase_density_n_lambda3=4*Ic/np.sqrt(np.pi),
          tidal_radius_over_R=rt,total_mass_Q=Qt))
    assert all(a['tidal_radius_over_R']<b['tidal_radius_over_R']
               for a,b in zip(rows[:-1],rows[1:]))
    assert all(a['total_mass_Q']<b['total_mass_Q']
               for a,b in zip(rows[:-1],rows[1:]))

    # Analytic monotonicity of Icritical, tested by finite differences.
    from critical_lowered import I
    for eta in (.5,1.,2.,4.,6.):
        eps=1e-4
        derivative=(I(eta+eps,eta+eps,nq=240)-I(eta-eps,eta-eps,nq=240))/(2*eps)
        exact=eta**1.5*np.exp(eta)/(3*np.expm1(eta)**2)
        assert abs(derivative/exact-1)<2e-5

    ideal_I=np.sqrt(np.pi)/4*zeta(1.5)
    ideal_slope=np.pi/(2*np.sqrt(6)*ideal_I)
    assert abs(ideal_slope-.553979545797)<1e-10
    assert rows[-1]['linear_center_slope']>ideal_slope
    assert abs(rows[4]['tidal_radius_over_R']-40.298792)<3e-5
    return dict(model="exactly saturated central lowered BE family psi0=eta",
      untruncated_Ic=ideal_I,untruncated_slope=ideal_slope,
      small_eta_Ic="(2/3)*sqrt(eta)-eta^(5/2)/90+...",
      small_eta_slope="3*pi/(4*sqrt(6*eta))",
      large_eta="Icritical tends to sqrt(pi)*zeta(3/2)/4; finite radius/mass diverge in limit",
      computed_profiles=rows,tests_passed=9,
      caveat="Mathematical cutoff family; no tidal kinetic model, galaxy fit or stability proof.")
if __name__=="__main__":print(json.dumps(run(),indent=2))
