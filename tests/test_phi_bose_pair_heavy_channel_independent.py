"""Independent swapped-mass Yukawa heavy χ1-antiχ1→2real φ validation."""
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from kinetics.phi_bose_pair_independent_vectorized import (
    sigma_independent,reaction_density_MB,ratio_to_H
)
from kinetics.phi_bose_pair_physical_onshell_quantum import differential_pair_sigma

def test_independent_swapped_dirac_mass_pair_CM_positive_spin_trace_and_sigma():
    z,w=leggauss(128)
    u=14*z;c=np.tanh(u);weight=14*w*(1-c*c)
    for E in (.21,.3,.5,1.,5.,100.):
        s=E*E
        independent=sigma_independent(s,128,m_in=.1,m_virtual=.02)
        from_lab_tree=float(np.dot(weight,differential_pair_sigma(
            s,c,m_in=.1,m_virtual=.02)))
        assert independent>0 and from_lab_tree>0
        assert math.isclose(independent,from_lab_tree,rel_tol=6e-8)
    assert sigma_independent(.1**2,128,m_in=.1,m_virtual=.02)==0.


def test_independent_swapped_heavy_pair_MB_thermal_integral_vs_previous_private_oracle():
    for T,target in ((.1,5.601313956908313e-21),
                     (1.,3.4912014442108976e-16)):
        e=reaction_density_MB(T,128,m_in=.1,m_virtual=.02)
        assert e["gamma_forward_GeV4"]>0
        assert math.isclose(e["gamma_forward_GeV4"],target,rel_tol=2e-5)
        assert e["mean_total_energy_per_event_GeV"]>=.2
    assert ratio_to_H(1.,61.75,128,m_in=.1,m_virtual=.02)>1000.


def test_heavy_pair_swapped_virtual_mass_remains_finite_under_angle_refinement():
    e=.5;s=e*e
    values=[sigma_independent(s,n,m_in=.1,m_virtual=.02)
            for n in (32,48,64,96,128)]
    assert all(x>0 for x in values)
    assert abs(values[-2]/values[-1]-1)<5e-7
    assert math.isclose(values[-1],
                        sigma_independent(s,160,m_in=.1,m_virtual=.02),
                        rel_tol=6e-8)
