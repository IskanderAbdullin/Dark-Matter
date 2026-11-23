"""Exact analytic energy integration of on-shell thermal FD/BE source."""
import math
from scipy.integrate import quad
from scipy.special import expit
from kinetics.phi_bose_thermal_exact_weak_reference import (
    daughter_bounds,scalar_thermal_number_integral,analytic_number_rate,
    gauss_weak,T0)


def test_analytic_FD_BE_daughter_integral_matches_independent_adaptive_energy():
    for p in (.1,1.,2.,8.,16.):
        E1,lo,hi,b=daughter_bounds(p)
        exact=scalar_thermal_number_integral(E1,lo,hi)
        test,_=quad(lambda k:expit(-(E1-k)/T0)/math.expm1(k/T0),
                    lo,hi,epsabs=1e-15,epsrel=1e-12)
        assert abs(exact/test-1)<2e-10


def test_analytic_1D_parent_number_rate_independent_of_positive_2D_gauss():
    for p in (2.,8.,16.):
        expected=analytic_number_rate(p)['Rphi']
        out=gauss_weak(p,96,48)
        assert out['Rphi']<0 and out['H']>0
        assert out['min_positive_event_weight']>0
        assert out['mean_Ephi_signed']>0
        assert abs(out['total_energy_defect'])<3e-25
        assert abs(out['Rphi']/expected-1)<2e-4


def test_three_resolution_sources_have_valid_energy_and_H_ledgers():
    for p in (2.,8.,16.):
        for nr,na in ((24,12),(48,24),(96,48)):
            a=gauss_weak(p,nr,na)
            assert a['Rphi']<0 and a['H']>0
            assert a['Qphi']<0 and a['Q1']>0 and a['Q2']<0
            assert abs(a['total_energy_defect'])<2e-24
