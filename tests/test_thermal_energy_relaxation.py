"""Independent FD thermodynamics and 5D collision thermal-response tests."""
import math
import sys
from pathlib import Path
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from thermal_energy_relaxation import (
    FD_thermo_moments,matched_mu,thermal_relaxation_rate)


def test_fixed_N_heat_capacity_and_mu_slope_vs_independent_fd_shift():
    M=.1;T=1.;h=1e-3
    ref=FD_thermo_moments(M,T,0)
    plus=matched_mu(M,T+h,ref["N"])
    minus=matched_mu(M,T-h,ref["N"])
    mu_slope=(plus-minus)/(2*h)
    assert math.isclose(mu_slope,ref["mu_slope_at_fixed_N"],
                        rel_tol=2e-6)
    rho_plus=FD_thermo_moments(M,T+h,plus)["rho"]
    rho_minus=FD_thermo_moments(M,T-h,minus)["rho"]
    capacity=(rho_plus-rho_minus)/(2*h)
    assert math.isclose(capacity,ref["c_N"],rel_tol=2e-6)
    assert capacity>0


def test_positive_cooling_and_dt_grid_refinement_not_spectral_gamma():
    a=thermal_relaxation_rate(delta_T=.04,n=16)
    b=thermal_relaxation_rate(delta_T=.01,n=24)
    assert a["Q_chi_minus_plus_GeV5"][0]>0
    assert a["Q_chi_minus_plus_GeV5"][1]<0
    assert a["gamma_T_energy_moment_GeV"]>0
    assert abs(b["gamma_T_energy_moment_GeV"]/
               a["gamma_T_energy_moment_GeV"]-1)<.002
    assert all(abs(x)<1e-10 for x in b["number_residual_minus_plus_fraction"])
    assert not b["is_true_spectral_eigenvalue_or_BGK_gamma"]


def test_same_T_neutrino_completion_conditional_scales_1p5():
    a=thermal_relaxation_rate(delta_T=.02,n=12,ge=4.)
    b=thermal_relaxation_rate(delta_T=.02,n=12,ge=6.)
    assert math.isclose(b["gamma_T_energy_moment_GeV"]/
                        a["gamma_T_energy_moment_GeV"],1.5,rel_tol=1e-12)


@pytest.mark.parametrize("parameters",[
    {"delta_T":.3},{"Te":0.},{"G":-1.},{"M":0.}
])
def test_reject_unphysical_input(parameters):
    with pytest.raises(ValueError):
        thermal_relaxation_rate(**parameters)
