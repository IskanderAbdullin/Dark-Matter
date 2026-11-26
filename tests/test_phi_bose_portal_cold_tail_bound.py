import math
from scipy.integrate import quad
from kinetics.phi_bose_portal_cold_tail_bound import row,ALPHA

def test_exact_tail_integral():
    for x in (2.,10.,50.):
        a=quad(lambda z:math.exp(-z)*(1+2/(x+z)+2/(x+z)**2),0,100)[0]
        assert math.isclose(a,1+2/x,rel_tol=1e-10)

def test_low_reheat_cold_count():
    assert .0096<ALPHA<.0097
    r=row(.02)
    assert 51<r["x_pair"]<53
    assert 2e-14<10**r["log10_omega_cold_upper"]<3e-14
    assert 10**row(.03)["log10_omega_cold_upper"]<.06
