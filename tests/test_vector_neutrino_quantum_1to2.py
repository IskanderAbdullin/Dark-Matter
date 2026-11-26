"""Fully quantum Zprime ↔ thermal ννbar 1<->2 exact angular controls."""
import math
import pytest
from kinetics.vector_neutrino_quantum_1to2 import (
 coeff_gauss,coeff_closed,mode_source,integrated_moments,ME_GEV)

@pytest.mark.parametrize("x,q",[(.03,0),(.03,.1),(.03,1),(.1,2),(.3,.05),(.3,5)])
def test_independent_closed_fd_angular_blocking(x,q):
    a=coeff_closed(x,q,1)
    b=coeff_gauss(x,q,1,nangle=96)
    assert math.isclose(a,b,rel_tol=2e-12,abs_tol=1e-13)
    assert 0<a<=1

def test_zero_momentum_FD_blocking_analytic_tanh():
    for x in (.03,.1,.3,1):
        assert math.isclose(coeff_closed(x,0,1),math.tanh(x/4),rel_tol=1e-15)

def test_mode_detailed_balance_and_inverse_production():
    m=.03;T=1;p=.6;g=.001
    eq=mode_source(m,p,T,g,1/math.expm1(math.hypot(p,m)/T))
    assert abs(eq["dbdt"])<1e-24
    for factor,sgn in [(.5,1),(2.,-1),(0.,1)]:
        a=mode_source(m,p,T,g,factor*eq["beq"])
        assert sgn*a["dbdt"]>0
        assert a["nu_energy_source_per_mode"]==-a["E"]*a["dbdt"]

def test_integrated_quantum_number_energy_and_entropy():
    m=.03*ME_GEV;T=ME_GEV;g=3e-7
    for factor,sgn in [(0,1),(.5,1),(2.,-1)]:
        a=integrated_moments(m,T,g,occupancy_factor=factor,n=128)
        assert sgn*a["dN_Zprime_GeV4"]>0
        assert a["dQ_Zprime_GeV5"]==-a["dQ_nu_GeV5"]
        assert a["dN_nu_plus_antinu_GeV4"]==-2*a["dN_Zprime_GeV4"]
        if factor>0:
            assert a["entropy_total_GeV4"]>=0
            assert math.isclose(a["entropy_Zprime_GeV4"]+a["entropy_nu_reservoir_GeV4"],
                                a["entropy_total_GeV4"],rel_tol=1e-10)
    eq=integrated_moments(m,T,g,occupancy_factor=1,n=128)
    assert abs(eq["dN_Zprime_GeV4"])<1e-30

def test_inverse_production_radial_convergence():
    m=.03*ME_GEV;T=ME_GEV;g=3e-7
    a=integrated_moments(m,T,g,occupancy_factor=0,n=128)
    b=integrated_moments(m,T,g,occupancy_factor=0,n=256)
    assert math.isclose(a["inverse_fill_rate_over_initial_equilibrium_density_GeV"],
                        b["inverse_fill_rate_over_initial_equilibrium_density_GeV"],rel_tol=1e-6)
