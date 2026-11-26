import math
from scipy.integrate import quad
from kinetics.phi_bose_portal_initial_gain_gate import integrated_s, M2, A, scalar_kinematics, M1, MCHI2, gain_initial

def test_resonance_integral():
    for lo,hi in [(1,30),(80,120),(200,500)]:
        m=.02
        actual=quad(lambda s:s*(s+2*m*m)/((s-M2)**2+A*A),lo,hi,points=[M2] if lo<M2<hi else None,epsrel=1e-10)[0]
        assert math.isclose(integrated_s(lo,hi,m),actual,rel_tol=1e-8)

def test_kinematics():
    e,lo,hi=scalar_kinematics(.012)
    assert hi>lo>=M1 and lo-e>=MCHI2

def test_reduced_gain():
    r=gain_initial(.6882487043,.012)
    assert abs(r["ratio"]-1.0002046815984291)<1e-6
    assert r["sharpened_gain_upper_reduced"]<2e-12
    assert r["sharpened_gain_upper_reduced"]<=r["integrated_positive_gain_upper_reduced"]
