"""Negative control: zero sampled cell source need not mean zero physical source."""
import math
import pytest
from kinetics.phi_bose_sparse_angle_false_zero import (
    audit, sampled_minimum, first_cell_edge,
)
from kinetics.phi_bose_ir_accessibility import min_over_parent_cutoff


@pytest.mark.parametrize("p,true_k,sample_k",[
    (2.,.00119925093757947,.009368256707822618),
    (8.,.000299988288329003,.009549351214434322),
    (16.,.000149998547511642,.009904824463004558),
])
def test_12_by_12_sparse_gauss_misses_physical_ir_support(p,true_k,sample_k):
    assert math.isclose(min_over_parent_cutoff(p),true_k,rel_tol=2e-9)
    assert math.isclose(sampled_minimum(p),sample_k,rel_tol=1e-11)
    assert sampled_minimum(p)>min_over_parent_cutoff(p)


def test_previous_n16_p2_first_cell_false_zero_not_kinematic():
    e=first_cell_edge(2.,16)
    assert math.isclose(e,.00860359519316437,rel_tol=1e-12)
    assert min_over_parent_cutoff(2.)<e<sampled_minimum(2.)
    row=next(x for x in audit() if x["pmax_GeV"]==2. and x["n"]==16)
    assert row["quadrature_false_zero"]
    assert row["physical_support_in_first_cell"]
    assert not row["sampled_event_in_first_cell"]


def test_n64_p2_first_cell_genuinely_kinematically_empty():
    row=next(x for x in audit() if x["pmax_GeV"]==2. and x["n"]==64)
    assert not row["physical_support_in_first_cell"]
    assert not row["quadrature_false_zero"]


def test_three_cutoffs_and_four_scalar_meshes():
    rows=audit()
    assert len(rows)==12
    assert sum(x["quadrature_false_zero"] for x in rows)>=2
    assert all(x["sampled_min_GeV"]>=x["physical_min_GeV"] for x in rows)
