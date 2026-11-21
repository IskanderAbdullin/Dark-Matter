"""Quantify nonuniform approach of bounded k=2 equilibria to critical Bose center.
Requires numpy/scipy and reproduce_cutoff_ladder.py; no GitHub Actions.
"""
import json, numpy as np
from math import pi,sqrt
from reproduce_cutoff_ladder import halo,I,U,W

def Iprime(psi,eta=4):
    q=sqrt(psi)*U
    x=eta-psi+q*q
    b=1/np.expm1(x)
    a=1/np.expm1(eta)*(1+1/np.expm1(eta))
    return psi**1.5*np.sum(W*U**2*(b*(1+b)-a))

def main():
    eta=4;Ic=float(I(eta,eta,2));theory=pi/(24*Ic)
    rows=[]
    for d in (.2,.1,.03,.01,.003,.001,.0003,.0001,3e-5,1e-5):
        p=eta-d;h=halo(2,p)
        curvature=Iprime(p)/(6*float(I(p,eta,2)))
        rows.append({'delta':d,
          'fmax':1/np.expm1(d)-1/np.expm1(eta)-p*np.exp(eta)/np.expm1(eta)**2,
          'radius_over_R':h['radius_over_R_or_limit'],
          'curvature_sqrt_delta':curvature*sqrt(d)})
    critical=halo(2,eta)
    assert abs(rows[-1]['curvature_sqrt_delta']/theory-1)<.03
    assert abs(rows[-1]['radius_over_R']/critical['radius_over_R_or_limit']-1)<.02
    assert rows[-1]['fmax']>90000
    return {'analytic_curvature_sqrt_delta_limit':theory,
      'critical_radius_over_R':critical['radius_over_R_or_limit'],
      'subcritical_to_critical':rows,'tests_passed':3,
      'caveat':'Bounded DF and its sup norm do not converge to the singular exact-critical endpoint.'}
if __name__=='__main__':print(json.dumps(main(),indent=2))
