import math
from kinetics.phi_bose_portal_recycling_energy_bound import KD,KINV,JCOEF,direct_sigma_v,energy_source
from kinetics.phi_bose_portal_low_reheat_gate import gamma

def test_conversion_bound_and_energy_moment():
    for s in (.040001,.05,.1,1.,100.,1e6):
        assert 0<=direct_sigma_v(s)<=KD
    T=.02;h=T*1e-5
    d=T*T*(gamma(T+h,.02)-gamma(T-h,.02))/(2*h)
    assert abs(energy_source(T,.02)/d-1)<1e-6
    assert KINV>2*KD
    assert .0005<JCOEF*.02**3<.00052
