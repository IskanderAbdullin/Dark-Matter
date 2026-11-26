"""Conservative heavy-drag matching and strict relativistic fail-closed."""
import math
import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from microphysics_nr_drag import heavy_nr_fokker_planck_coefficients


def test_einstein_relation_and_independent_heavy_asymptote():
    M=100.;T=.001;G=.2
    r=heavy_nr_fokker_planck_coefficients(M,T,G)
    Q=r["Q_rest_GeV2_per_chi"]
    assert Q>0
    assert math.isclose(r["leading_drag_eta_GeV"],Q/(3*T),rel_tol=1e-13)
    assert math.isclose(r["momentum_diffusion_kappa_GeV3"]/
                        (2*M*r["leading_drag_eta_GeV"]),T,rel_tol=1e-13)
    assert abs(r["heavy_limit_relative_deviation"])<.001


def test_yukawa_free_elastic_zero_portal_zero_drag():
    r=heavy_nr_fokker_planck_coefficients(100.,.001,0.)
    assert r["leading_drag_eta_GeV"]==0
    assert r["momentum_diffusion_kappa_GeV3"]==0


def test_relativistic_parent_must_not_get_fake_BGK_gamma():
    with pytest.raises(ValueError,match="Relativistic"):
        heavy_nr_fokker_planck_coefficients(.1,1.,.5)
    with pytest.raises(ValueError):
        heavy_nr_fokker_planck_coefficients(1.,0.,.5)
    with pytest.raises(ValueError):
        heavy_nr_fokker_planck_coefficients(100.,.001,.2,min_mass_over_temperature=10.)
