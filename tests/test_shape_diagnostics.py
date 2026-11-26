"""Independent cumulative mass/quantile checks for pre-freezeout Bose spectra."""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from shape_diagnostics import (
    validate_snapshot,cdf_number_at_q,quantile_number_q,
    describe_number_spectrum,one_case,run
)


def uniform_bose_spectrum():
    edges=np.array([0.,1.,2.,3.])
    number=(edges[1:]**3-edges[:-1]**3)/(6*np.pi**2)
    return {"q_bin_edges":edges.tolist(),
            "phi_occupations":[1.,1.,1.],
            "phi_number_per_bin":number.tolist(),
            "total_phi_number":float(number.sum())}


def test_exact_q2_phase_volume_cdf_of_uniform_population():
    data=uniform_bose_spectrum()
    edges,b,n,total=validate_snapshot(data)
    assert cdf_number_at_q(edges,n,0)==0
    assert cdf_number_at_q(edges,n,10)==1
    assert abs(cdf_number_at_q(edges,n,1.5)-.125)<1e-12
    assert abs(quantile_number_q(edges,n,.5)-3*(.5)**(1/3))<1e-12
    assert abs(quantile_number_q(edges,n,.1)-3*(.1)**(1/3))<1e-12


def test_all_probabilities_produce_monotone_quantiles():
    data=uniform_bose_spectrum()
    e,_,n,_=validate_snapshot(data)
    vals=[quantile_number_q(e,n,p) for p in (0,.01,.1,.5,.9,.99,1)]
    assert vals==sorted(vals)
    assert vals[0]==0.
    assert vals[-1]==3.


def test_empty_and_invalid_spectrum_rejected():
    data=uniform_bose_spectrum()
    data["phi_number_per_bin"]=[0,0,0]
    data["total_phi_number"]=0
    with pytest.raises(ValueError):
        describe_number_spectrum(data)
    with pytest.raises(ValueError):
        quantile_number_q([0,1],[1],1.1)


def test_short_closed_schur_snapshot_conserves_bin_number():
    case=one_case(5000,tend=.0002)
    d=case["summary"]
    assert d["bin_count"]==199
    assert d["number_over_T3"]>0
    assert 0<=d["fraction_q_below_0p1"]<=1
    qs=d["quantile_momenta"]
    assert 0<qs["q01"]<=qs["q10"]<=qs["q50"]<=qs["q90"]<=qs["q99"]
    values=list(d["normalized_number_CDF"].values())
    assert values==sorted(values)
    assert max(abs(v) for v in case["number_and_energy_residuals"].values())<1e-8
    assert case["short_time_only"]


def test_quick_shape_survey_explicitly_refuses_class_status():
    result=run(quick=True)
    assert result["shape_difference"] is None
    assert "NOT COSMOLOGICAL" in result["status"]
    assert "do not pass to CLASS" in result["prohibition"]
