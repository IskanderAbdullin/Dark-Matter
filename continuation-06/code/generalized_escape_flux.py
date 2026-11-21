"""Generalized near-escape power-law diffusion in a toy 1D Bose flux model.
J=D(s)[df/ds-beta*f*(1+f)], D(s)=D0*s**alpha.
NOT a derived orbit-averaged C22 or physical galactic evaporation rate.
"""
from __future__ import annotations
import json
import numpy as np
from scipy.integrate import solve_ivp

def solve_one(alpha,beta=1.0,J=0.01,D0=1.0,smin=1e-8,smax=0.04):
    assert alpha<1
    k=1-alpha
    f0=J/(D0*k)*smin**k
    sol=solve_ivp(lambda s,f:[J/D0*s**(-alpha)+beta*f[0]*(1+f[0])],
        [smin,smax],[f0],rtol=1e-11,atol=1e-23,
        max_step=.0002,dense_output=True)
    if not sol.success:raise RuntimeError(sol.message)
    s=np.array([1e-6,1e-5,1e-4,.001,.01,.03])
    f=sol.sol(s)[0]
    asym=J/(D0*k)*s**k
    ratios=f/asym
    assert np.all(f>0) and np.all(np.diff(f)>0)
    assert abs(ratios[0]-1)<.005
    return dict(alpha=float(alpha),
      predicted_edge_power=float(k),
      numerical_log_slope=float(np.polyfit(np.log(s[:4]),np.log(f[:4]),1)[0]),
      ratios_to_analytic_leading=list(map(float,ratios)),
      sampled_binding=list(map(float,s)),
      residual_at_smallest=float(ratios[0]-1))

def run():
    cases=[solve_one(a) for a in (-2,-1,-.5,0,.5,.8)]
    for x in cases:
        assert abs(x["numerical_log_slope"]-x["predicted_edge_power"])<.03
    return dict(model="1D toy Bose energy flux J=D(s)[f_s-beta*f(1+f)]",
      assumed_positive_escape_flux=.01,beta=1,D0=1,
      no_finite_absorbing_solution_if_alpha_ge1=True,
      scenarios=cases,checks_passed=2*len(cases),
      caveat="No real C22 collision kernel, energy-state diffusion, dynamical galaxy escape or observable rate.")

if __name__=="__main__":print(json.dumps(run(),indent=2))
