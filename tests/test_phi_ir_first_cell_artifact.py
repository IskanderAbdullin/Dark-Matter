"""First Bose energy bin is NOT a physical condensate: finite p²dp test."""
import math
import pytest
from scipy.special import zeta
from kinetics.phi_ir_first_cell_artifact import (
    thermal_number_integral,thermal_first_cell_and_total,
)


def test_massless_thermal_BE_has_finite_integrated_number_without_condensate():
    exact=zeta(3)/math.pi**2
    got=thermal_number_integral(0,1,50)
    assert abs(got/exact-1)<1e-12
    with pytest.raises(ValueError):
        thermal_first_cell_and_total(0,bins=64)


def test_toy_mphi_overestimates_integrated_first_bose_cell():
    x=thermal_first_cell_and_total(.001,bins=64)
    assert 60<x['first_cell_false_enhancement_ratio']<70
    assert 1.1<x['full_false_enhancement_ratio']<1.2
    assert abs(x['only_first_cell_corrected_ratio']-1)<.004


def test_86eV_example_is_a_finite_cell_artifact_not_a_BEC():
    x=thermal_first_cell_and_total(86e-9,bins=64,T=1.)
    assert x['first_cell_false_enhancement_ratio']>7e5
    assert x['full_false_enhancement_ratio']>1e3
    assert 0<x['first_cell_exact_thermal_number']<.001
    assert abs(x['only_first_cell_corrected_ratio']-1)<.004
    assert abs(x['first_cell_false_enhancement_ratio']/x['asymptotic_first_cell_ratio_h_over_3m']-1)<.04


def test_first_cell_correction_improves_total_at_three_meshes_without_thermostat():
    e=[]
    for n in (16,64,256):
        x=thermal_first_cell_and_total(86e-9,bins=n)
        assert x['full_false_enhancement_ratio']>10
        e.append(abs(x['only_first_cell_corrected_ratio']-1))
    assert e[0]>e[1]>e[2]
    assert e[-1]<.0004


def test_negative_nonphysical_thermal_integral_inputs_rejected():
    for m,T,p in ((-1,1,3),(0,0,3),(0,1,-1)):
        with pytest.raises(ValueError):
            thermal_number_integral(m,T,p)
