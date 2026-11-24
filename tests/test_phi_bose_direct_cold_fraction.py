"""Independent physical on-shell unenhanced Yukawa direct cold fraction tests.

Not a late-time cosmological relic; direct FD χ1 decays only, no f2,bφ,
inverse, plasma or Friedmann. Ultra-relativistic closed form is independent
from numerical exact massive event/angle integration.
"""
import math

import numpy as np
from scipy.integrate import quad
from scipy.special import zeta

from kinetics.phi_bose_direct_cold_fraction import (
    ultra_rel_kinematic_fraction,ultra_rel_number_pdf,
    ultra_rel_mean_k_over_T,finite_T_direct_cold_fraction,
    finite_T_direct_mean_x,
)


def test_ultrarelativistic_analytic_pdf_cdf_and_mean_independently_integrated():
    A=.96
    norm=quad(ultra_rel_number_pdf,0.,100.,epsabs=1e-13)[0]
    first_moment=quad(lambda x:x*ultra_rel_number_pdf(x),
                      0.,100.,epsabs=1e-13)[0]
    assert math.isclose(norm,1.,rel_tol=1e-12)
    assert math.isclose(first_moment,ultra_rel_mean_k_over_T(),rel_tol=1e-12)
    assert math.isclose(first_moment,9*A*zeta(3.,1.)/math.pi**2,rel_tol=1e-12)
    for threshold in (.001,.01,.1,1.,3.):
        cdf=quad(ultra_rel_number_pdf,0,threshold,
                 epsabs=1e-13)[0]
        assert math.isclose(
            cdf,ultra_rel_kinematic_fraction(threshold),
            abs_tol=2e-13,rel_tol=2e-11)
    assert math.isclose(ultra_rel_kinematic_fraction(.01),
                        .008745886700982197,rel_tol=5e-12)


def test_massive_positive_exact_angle_production_converges_to_ultrarel_limit():
    for alpha in (.001,.01,.1):
        reference=ultra_rel_kinematic_fraction(alpha)
        out=finite_T_direct_cold_fraction(1000.,alpha)
        assert out["cold_direct_decay_fraction"]>0
        assert abs(out["cold_direct_decay_fraction"]/reference-1)<3e-5
        assert out["physical_parent_p_min_GeV"]>0
    assert math.isclose(finite_T_direct_mean_x(1000.),
                        1.0522986971538855,rel_tol=5e-9)
    assert math.isclose(ultra_rel_mean_k_over_T(),
                        1.0522986759380715,rel_tol=5e-12)


def test_parent_boost_and_temperature_change_cold_fraction_enormously():
    # For fixed massive Dirac model and x=k/T<.01, ordinary unenhanced
    # direct decays at 0.1GeV rely on a rare parent tail (p>=2.4GeV).
    r01=finite_T_direct_cold_fraction(.1,.01)
    r1=finite_T_direct_cold_fraction(1.,.01)
    r100=finite_T_direct_cold_fraction(100.,.01)
    assert math.isclose(r01["cold_direct_decay_fraction"],
                        2.7658668662688188e-14,rel_tol=2e-5)
    assert math.isclose(r1["cold_direct_decay_fraction"],
                        .00521753296596788,rel_tol=3e-7)
    assert math.isclose(r100["cold_direct_decay_fraction"],
                        .008744159100937229,rel_tol=3e-7)
    assert 2.39<r01["physical_parent_p_min_GeV"]<2.41
    assert r1["cold_direct_decay_fraction"]/r01["cold_direct_decay_fraction"]>1e11
    assert r100["cold_direct_decay_fraction"]<ultra_rel_kinematic_fraction(.01)


def test_tiny_nonzero_thermal_tail_not_falsely_declared_kinematic_zero():
    # If the parent grid is hard-truncated at q=p/T=75, T=.1 and
    # x=k/T<.001 would look exactly zero, but exact on-shell parent
    # threshold is q~240 and the FD tail is strictly nonzero.
    x=finite_T_direct_cold_fraction(.1,.001)
    assert 239.<x["required_parent_q_min"]<241.
    assert 0<x["cold_direct_decay_fraction"]<1e-100
    assert math.isclose(x["cold_direct_decay_fraction"],
                        4.6567561619051507e-110,rel_tol=3e-4)
    assert x["parent_q_max_for_numerator"]>300.


def test_unenhanced_is_not_self_consistent_Bose_burst_relic():
    # Source is forward only: no inverse, stimulation, SM or FRW.
    # Even in the strict high-T limit less than 1% of initial direct
    # scalar daughters have k/T<.01.
    assert 0<ultra_rel_kinematic_fraction(.001)<.001
    assert 0<ultra_rel_kinematic_fraction(.01)<.01
    assert .08<ultra_rel_kinematic_fraction(.1)<.09
    assert ultra_rel_mean_k_over_T()>1.
    assert ultra_rel_mean_k_over_T()<2.701178032919
