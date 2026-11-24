"""Independent vectorized Yukawa pair σ and MB thermal collision checks."""
import math
from kinetics.phi_bose_pair_annihilation_mb import (
    cross_section,thermal_forward_rate_MB
)
from kinetics.phi_bose_pair_mb_detailed_balance import (
    equilibrium_pair_densities
)
from kinetics.phi_bose_pair_independent_vectorized import (
    sigma_independent,reaction_density_MB,ratio_to_H
)

def test_independent_vectorized_minkowski_trace_and_phi_half_factor():
    for sqrt_s in (.05,.1,.5,1.,5.,100.):
        s=sqrt_s**2
        slow=cross_section(s,nangle=128)
        fast=sigma_independent(s,angle_order=128)
        assert fast>0
        assert math.isclose(fast,slow,rel_tol=5e-8)

def test_independent_angular_refinement_is_needed_for_thermal_pair_rate():
    for T in (.1,1.,100.):
        values=[reaction_density_MB(T,n)["gamma_forward_GeV4"]
                for n in (32,48,64,96,128)]
        assert all(v>0 for v in values)
        assert abs(values[-2]/values[-1]-1)<1e-6
        # Not all 48 point angular integrations reach <1e-4 at T=.1.
        if T==.1:
            assert abs(values[1]/values[-1]-1)>7e-4
        if T==1.:
            assert abs(values[1]/values[-1]-1)>3e-4

def test_gamma_T_and_H_independent_MBoltzmann_reference():
    for T,g,expected in ((.1,17.25,6936.017005351313),
                          (1.,61.75,964.7820572539009),
                          (100.,106.75,23.40225556833799)):
        fast=reaction_density_MB(T,128)
        slow=thermal_forward_rate_MB(T,g,nangle=128)
        assert math.isclose(fast["gamma_forward_GeV4"],
                            slow["forward_pair_reaction_density_GeV4"],rel_tol=1e-5)
        assert math.isclose(ratio_to_H(T,g),expected,rel_tol=1e-6)

def test_positive_pair_energy_moment_matches_beta_finite_difference():
    T=1.
    e=reaction_density_MB(T,128)
    d=.003
    plus=reaction_density_MB(T+d,128)["gamma_forward_GeV4"]
    minus=reaction_density_MB(T-d,128)["gamma_forward_GeV4"]
    q_num=T*T*(plus-minus)/(2*d)
    assert e["pair_energy_forward_GeV5"]>0
    assert math.isclose(q_num,e["pair_energy_forward_GeV5"],rel_tol=2e-4)
    assert math.isclose(e["mean_total_energy_per_event_GeV"],
                        4.409061907291514,rel_tol=5e-6)
    eq=equilibrium_pair_densities(T,128)
    assert math.isclose(eq["gamma_forward_equilibrium_GeV4"],
                        e["gamma_forward_GeV4"],rel_tol=2e-5)
