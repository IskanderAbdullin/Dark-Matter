"""Independent continuum rectangular phase-space cutoff tests."""
from __future__ import annotations
import math
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import dblquad
from scipy.special import expit

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from uv_source_decomposition import exact_rectangle_source,one_comparison,run
from onshell_collision_2d import vacuum_reference


@pytest.mark.parametrize("Q",[.1,.2,.5,1.,3.])
def test_rectangle_analytic_daughter_integral_matches_independent_2d_quadrature(Q):
    m=.1
    c=m*m/4.
    pref=4*m/(2*math.pi**2)
    direct=dblquad(lambda chi,phi:float(expit(-(phi+chi))),
                   c/Q,Q,lambda phi:c/phi,lambda phi:Q,
                   epsabs=1e-12,epsrel=1e-10)[0]*pref
    analytic=exact_rectangle_source(Q,mass=m)
    assert math.isclose(direct,analytic,rel_tol=2e-10,abs_tol=1e-13)


def test_no_kinematic_pairs_below_threshold():
    assert exact_rectangle_source(.049)==0.
    assert exact_rectangle_source(.05)==0.
    assert exact_rectangle_source(.1)>0
    with pytest.raises(ValueError):
        exact_rectangle_source(0)


def test_rectangle_rate_increases_to_full_continuum_limit():
    results=[exact_rectangle_source(q) for q in (.5,2.,5.,7.,10.,20.,30.)]
    assert all(a<b for a,b in zip(results,results[1:]))
    assert abs(results[-1]/vacuum_reference()-1)<1e-10


@pytest.mark.parametrize("maximum",[5000,7000,8500,10000])
def test_error_decomposition_identity_and_continuum_domain(maximum):
    d=one_comparison(maximum)
    assert abs(d["exact_multiplicative_decomposition_residual"])<1e-12
    assert d["continuum_infinite_domain_source"]>=d["continuum_finite_rectangle_source"]
    assert d["daughter_energy_max_EDGE"]>maximum*.001
    assert d["relative_cell_rate_error_within_same_domain"]<0
    assert not d["not_a_nonlinear_source_convergence_test"] is False


def test_cutoff_dominates_coarse_case_but_cell_error_plateaus_high_Q():
    a=one_comparison(5000)
    b=one_comparison(8500)
    c=one_comparison(10000)
    assert abs(a["relative_domain_truncation_error"])>.01
    assert abs(b["relative_domain_truncation_error"])<.001
    assert abs(c["relative_domain_truncation_error"])<.0002
    assert abs(c["relative_cell_rate_error_within_same_domain"])>.00009
    assert abs(c["relative_cell_rate_error_within_same_domain"])<.00012
    assert abs(b["relative_cell_rate_error_within_same_domain"]-
               c["relative_cell_rate_error_within_same_domain"])<2e-6


def test_quick_results_are_labeled_vacuum_only():
    out=run(quick=True)
    assert len(out["cases"])==2
    assert out["alternative_meshing_at_same_approximate_Q"]==[]
    assert "NOT" not in out["main_result"]
    assert "nonlinear" in out["caveat"]
