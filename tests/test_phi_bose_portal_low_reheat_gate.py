"""Low reheating contact and finite-mass checks."""
from math import pi
from kinetics.phi_bose_portal_low_reheat_gate import gamma, COUP2, M, LEPTONS

def test_contact_and_tail():
    T=.001
    expected=COUP2*T**8*384/(48*pi**5*M**4)
    assert abs(gamma(T,0.,0.)/expected-1)<1e-5
    for T in (.005,.03,.1):
        for m in (.02,.1):
            a=gamma(T,m,0.,100)
            b=gamma(T,m,0.,150)
            assert abs(a/b-1)<1e-7
            assert gamma(T,m,LEPTONS[0])<=a*(1+1e-9)
