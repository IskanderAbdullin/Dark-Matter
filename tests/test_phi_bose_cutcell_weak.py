"""Regression for physically nonzero first Bose cell missed by 12x12 Gauss."""
import math
from kinetics.phi_bose_cutcell_weak import audit, integrate
from kinetics.phi_bose_sparse_angle_false_zero import first_cell_edge


def test_first_cell_exact_cut_angle_integral_is_nonzero():
    a=audit()
    assert math.isclose(a["first_cell_parent_threshold_GeV"],
                        .269991073952339,rel_tol=3e-9)
    assert math.isclose(a["first_cell_Rphi_GeV4"],
                        -2.762149360447387e-12,rel_tol=3e-6)
    assert math.isclose(a["first_cell_Qphi_GeV5"],
                        -1.2277997619515968e-14,rel_tol=3e-6)
    assert a["first_cell_Rphi_GeV4"]<0
    assert .20<a["first_cell_fraction_of_total"]<.23


def test_total_full_boosted_angle_12x12_underestimates_rate():
    a=audit()
    assert math.isclose(a["total_Rphi_GeV4"],
                        -1.2954539554916677e-11,rel_tol=3e-6)
    assert .09<a["old_12x12_fractional_total_underestimate"]<.11


def test_cutcell_is_strictly_zero_when_entire_cell_below_true_support():
    e=first_cell_edge(2.,64)
    assert integrate(2.,e)==0.
    assert integrate(2.,e,True)==0.
