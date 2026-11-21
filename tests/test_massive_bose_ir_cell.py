"""Deep-IR massive Bose occupation source is not a condensate number."""
from __future__ import annotations
import math
import sys
from pathlib import Path

import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from massive_bose_ir_cell import initial_ir_cell,run


def test_very_large_qzero_occupation_corresponds_to_tiny_phase_space_number():
    row=initial_ir_cell(M=.1,mphi=.001,q_cut=.0001)
    assert 26000<row["exact_qzero_occupation_rate"]<28000
    assert 26000<row["exact_cell_phase_volume_averaged_occupation_rate"]<28000
    assert 1.<row["cell_average_to_qzero_source_ratio"]<1.02
    assert 2e-8<row["fraction_all_initial_phi_decays_below_qcut"]<4e-8
    assert row["is_a_condensate_n0_source"] is False


@pytest.mark.parametrize("mass",[.0005,.001,.002,.005])
def test_cold_cell_mean_converges_to_exact_zero_momentum_occupation(mass):
    small=initial_ir_cell(mphi=mass,q_cut=1e-5)
    smaller=initial_ir_cell(mphi=mass,q_cut=1e-6)
    assert abs(smaller["cell_average_to_qzero_source_ratio"]-1)<.0001
    assert abs(smaller["cell_average_to_qzero_source_ratio"]-1)<abs(
        small["cell_average_to_qzero_source_ratio"]-1)


def test_phase_space_number_scales_as_q_cubed_at_fixed_positive_phi_mass():
    rows=[initial_ir_cell(mphi=.002,q_cut=cut)
          for cut in (1e-5,2e-5,4e-5)]
    ratios=[rows[k+1]["exact_physical_number_flux_below_qcut"]/
            rows[k]["exact_physical_number_flux_below_qcut"]
            for k in (0,1)]
    assert all(7.9<x<8.1 for x in ratios)


def test_initial_ir_volume_can_be_small_even_with_very_large_occupation():
    a=initial_ir_cell(mphi=.002,q_cut=.0001)
    b=initial_ir_cell(mphi=.005,q_cut=.0001)
    assert a["exact_qzero_occupation_rate"]>35000
    assert b["exact_qzero_occupation_rate"]>8000
    assert a["fraction_all_initial_phi_decays_below_qcut"]<1e-7
    assert b["fraction_all_initial_phi_decays_below_qcut"]<1e-7


@pytest.mark.parametrize("bad",[
    {"q_cut":0},{"q_cut":-1},{"gphi":0},{"mphi":-.01}
])
def test_invalid_q_cutoff_and_mass_are_rejected(bad):
    with pytest.raises(ValueError):
        initial_ir_cell(**bad)


def test_survey_is_not_labeled_as_cosmological_condensate():
    record=run()
    assert len(record["cases"])==12
    assert "NOT VALIDATED" in record["zero_mode_or_FRW_classification"]
    assert all(row["is_a_condensate_n0_source"] is False for row in record["cases"])
