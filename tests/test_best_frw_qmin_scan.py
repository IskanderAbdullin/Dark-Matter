"""FRW qmin scan is a *real* full 3-species BEST diagnostic, not condensate."""
from __future__ import annotations
import math
import sys
from pathlib import Path

import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from best_frw_qmin_scan import qmin_protocol,execute_case
from radiation_collisionless_parent_bose import evaluate


def test_real_frw_scan_changes_only_comoving_ir_cutoff():
    cases=[qmin_protocol(qmin=cut) for cut in (1e-4,5e-5,2e-5)]
    assert all(x["kinetic_source"]["grid_n"]==96 for x in cases)
    assert all(x["steps"]==128 and x["dtGamma"]==3.125e-7 for x in cases)
    assert all(x["end_tGamma"]==4e-5 for x in cases)
    assert all(math.isclose(x["rd_final_expected_a"],math.sqrt(1.4))
               for x in cases)
    assert [x["kinetic_source"]["qmin"] for x in cases]==[
        1e-4,5e-5,2e-5]
    assert all(x["qmin_scan_is_not_macroscopic_BEC"] for x in cases)
    assert not any(x["nonzero_BEC_condensate_inferred"] for x in cases)


def test_independent_redshifted_Bose_regular_cell_scales_q_cubed():
    ref=[evaluate(qmin=cut)["below_qmin_comoving_number"]
         for cut in (1e-4,5e-5,2e-5)]
    assert math.isclose(ref[0],1.7092584190101826e-14,rel_tol=2e-7)
    assert math.isclose(ref[0]/ref[1],8,rel_tol=.015)
    assert math.isclose(ref[1]/ref[2],15.625,rel_tol=.006)


@pytest.mark.parametrize("bad",[
    {"qmin":0},{"qmin":-1},{"grid_n":2},
    {"neval":100},{"steps":0},{"dt":0},{"t0":0}
])
def test_bad_qmin_or_time_or_grid_are_rejected(bad):
    with pytest.raises(ValueError):
        qmin_protocol(**bad)


def test_missing_thirdparty_best_must_not_forge_a_series(tmp_path):
    out=tmp_path/"false_data.json"
    with pytest.raises(ValueError,match="missing"):
        execute_case(tmp_path,out,qmin_protocol(steps=1))
    assert not out.exists()
