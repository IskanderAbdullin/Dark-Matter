"""Focus on stiff non-asymptotic vs asymptotic second-order dt regime."""
import pytest
from kinetics.vector_neutrino_strang_dt_convergence import scan


def test_stiff_h0p1_coarse_negative_and_refined_order_two():
    r=scan(h=.1,tau=1.,dts=(.125,.0625,.03125,.015625,
                             .0078125,.00390625,.001953125))
    diffs=r["step_halving_differences"]
    assert diffs[1]["successive_error_ratio"]<1.0
    assert diffs[2]["successive_error_ratio"]<2.0
    assert diffs[-1]["successive_error_ratio"]>3.7
    assert diffs[-2]["successive_error_ratio"]>3.5
    assert all(row["min_step_entropy_gain"]>0 for row in r["rows"])
    assert all(row["max_abs_charge_or_energy_residual"]<2e-12
               for row in r["rows"])
    assert all(row["min_b"]>=0 and row["min_fermion"]>=0
               and row["max_fermion"]<=1 for row in r["rows"])


@pytest.mark.parametrize("h",(.4,.2))
def test_less_stiff_graph_enters_second_order_temporal_regime(h):
    r=scan(h=h,tau=1.,dts=(.0625,.03125,.015625,.0078125))
    ratios=[x["successive_error_ratio"] for x in
            r["step_halving_differences"][1:]]
    assert min(ratios)>3.0
