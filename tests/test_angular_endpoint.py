"""Angular endpoint source-support and false IR-regulator diagnostics."""
from __future__ import annotations
import sys
from pathlib import Path
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from angular_endpoint_audit import analytic_angular_floor, quadrature_support


@pytest.mark.parametrize("nmu",[48,72,96,120,200])
def test_analytic_floor_is_below_real_quadrature_support(nmu):
    d=quadrature_support(nmu,np1=140,nq=100)
    assert d["analytic_angular_support_floor_q"]<=d["lowest_actual_quadrature_daughter_q"]*(1+1e-10)
    assert d["qmin_below_all_generated_daughter_events"]
    assert d["angular_support_floor_over_nominal_qmin"]>50
    assert d["analytic_angular_support_floor_q"]>d["true_continuum_angle_floor_for_p_le_26"]


def test_angular_floor_decreases_with_angle_refinement():
    floor=[analytic_angular_floor(n)["analytic_angular_support_floor_q"] for n in (48,72,96,120,200)]
    assert all(a>b>0 for a,b in zip(floor,floor[1:]))
    assert floor[-1]<floor[0]/3


def test_angular_endpoint_is_independent_of_momentum_bin_count():
    a=quadrature_support(96,np1=140,nq=50)
    b=quadrature_support(96,np1=140,nq=150)
    assert a["lowest_actual_quadrature_daughter_q"]==b["lowest_actual_quadrature_daughter_q"]


def test_angular_floor_is_influenced_by_parent_momentum_cap():
    result=analytic_angular_floor(800)
    assert result["parent_p_optimal_at_mu"]>26
    assert result["parent_p_used_after_pmax_bound"]==26.
    assert result["analytic_angular_support_floor_q"]>result["true_continuum_angle_floor_for_p_le_26"]
