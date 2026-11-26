import numpy as np
from kinetics.phi_bose_portal_fd_pole_gate import angular_factor,ratio,rate

def test_quantum_factor_and_nwa_normalization():
    a=angular_factor(5.1,np.array([5.1,10.,20.]),100.,.02,1.)
    assert np.all((a>=.25)&(a<=1))
    assert abs(rate(1.,.02,False)/3.447115846e-16-1)<1e-7

def test_inversion_and_quadrature():
    a=ratio(5.,.0001,True,32)
    b=ratio(5.,.0001,True,64)
    assert abs(a["ratio"]-b["ratio"])<1e-9
    assert abs(b["ratio"]-1.0000637212455268)<1e-7
