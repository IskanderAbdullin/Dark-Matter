"""Analytic non-equilibrium Bose/Pauli protocol and external-execution gating."""
from __future__ import annotations
import math
import sys
from pathlib import Path
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from best_chemical_affinity_reference import affinity_protocol,execute_best,BEST_PIN


def test_chemical_affinity_external_protocol_matches_original_2d_integrals():
    p=affinity_protocol()
    assert p["best_pinned_sha"]==BEST_PIN
    assert not p["external_best_executed"]
    assert not p["self_consistent_frw_or_class"]
    assert p["chemical_affinity"]>0
    assert math.isclose(p["exact_analytic_number_rate"],.0008565027253502168,rel_tol=1e-8)
    E=p["exact_analytic_three_energy_rates"]
    assert math.isclose(E["chi1"],-.0018101982860539644,rel_tol=1e-7)
    assert math.isclose(E["phi"],.0006730451878862239,rel_tol=1e-7)
    assert math.isclose(E["chi2"],.0011371530981677406,rel_tol=1e-7)
    assert abs(sum(E.values()))<1e-8
    assert abs(p["independent_2D_energy_check"])<1e-8
    assert math.isclose(p["exact_analytic_entropy_production"],
                        .05*p["exact_analytic_number_rate"],rel_tol=1e-12)


def test_inverse_dominated_protocol_has_negative_net_source_and_positive_entropy():
    p=affinity_protocol(mu1=-.25)
    assert p["chemical_affinity"]<0
    assert p["exact_analytic_number_rate"]<0
    assert p["exact_analytic_entropy_production"]>0
    for name,rate in p["exact_analytic_three_energy_rates"].items():
        assert rate>0 if name=="chi1" else rate<0


@pytest.mark.parametrize("bad",[
    {"mu_phi":0.},{"mu_phi":.2},{"mu1":10.},
    {"grid_n":3},{"vegas_neval":100},{"mu1":float("nan")}
])
def test_unphysical_or_unreproducible_protocol_is_rejected(bad):
    with pytest.raises(ValueError):
        affinity_protocol(**bad)


def test_external_run_does_not_fake_missing_best_source(tmp_path):
    with pytest.raises(ValueError,match="not found"):
        execute_best(tmp_path,tmp_path/"fake.json")
    assert not (tmp_path/"fake.json").exists()
