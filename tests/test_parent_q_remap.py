"""Conservative moving-grid Liouville remap tests."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from parent_q_remap import rebin,number,energy,entropy,moving_parent_edges,thermal_shell_averages,diagnose,run


def test_exact_identity_constant_and_analytic_one_bin_average():
    old=[0.,1.,2.]
    assert np.array_equal(rebin(old,old,[.1,.6]),[.1,.6])
    assert np.isclose(rebin(old,[0.,2.],[1.,0.])[0],1/8)
    assert np.allclose(rebin(old,[0.,.3,1.5,2.],[.4,.4]),.4)


def test_random_remaps_preserve_number_pauli_and_increase_entropy():
    rng=np.random.default_rng(2026)
    for n1,n2 in ((6,20),(50,100),(100,50),(300,300)):
        a=np.r_[0.,np.sort(rng.uniform(.001,5.,n1-1)),5.]
        b=np.r_[0.,np.sort(rng.uniform(.001,5.,n2-1)),5.]
        f=rng.uniform(0,1,n1)
        mapped=rebin(a,b,f)
        assert 0<=mapped.min()<=mapped.max()<=1
        assert abs(number(a,f)-number(b,mapped))<1e-11
        assert entropy(b,mapped)>=entropy(a,f)-1e-10


def test_moving_parent_momentum_edges_fail_if_mass_crosses_reference_center():
    E=[.1,.13,.2,.4]
    q0=moving_parent_edges(E,.09,1.)
    q1=moving_parent_edges(E,.095,1.)
    assert q0[0]==q1[0]==0 and q0[-1]==q1[-1]==1
    assert np.all(q1[1:-1]<q0[1:-1])
    with pytest.raises(ValueError,match="rebuild resonant graph"):
        moving_parent_edges(E,.12,1.)


def test_thermal_shell_average_has_physical_fermi_bounds():
    edges=[0.,.01,.03,.1,.2,.5,1.,3.]
    f=thermal_shell_averages(edges,.1)
    assert np.all((f>0)&(f<1))
    report=diagnose(edges,[0.,.02,.04,.09,.3,.6,3.],f,.12)
    assert abs(report["absolute_number_error"])<1e-11
    assert report["min_fermion_occupation"]>=0
    assert report["max_fermion_occupation"]<=1
    assert report["numerical_fermion_entropy_increase"]>=0
    assert np.isfinite(report["relative_energy_regrid_error"])


def test_energy_projection_error_is_reported_instead_of_hidden():
    old=[0.,.5,1.,2.,4.]
    new=[0.,.3,1.2,1.8,4.]
    f=[1.,0.,.8,0.]
    report=diagnose(old,new,f,.5)
    assert abs(report["relative_energy_regrid_error"])>1e-5
    assert report["numerical_fermion_entropy_increase"]>0


def test_invalid_edges_occupations_or_domain_rejected():
    with pytest.raises(ValueError):rebin([0,1],[0,.8],[.2])
    with pytest.raises(ValueError):rebin([0,1],[0,1],[1.1])
    with pytest.raises(ValueError):rebin([0,1,1],[0,1],[.1,.1])
    with pytest.raises(ValueError):thermal_shell_averages([0,1],0)
    with pytest.raises(ValueError):energy([0,1],[.3],0)


def test_actual_resonant_parent_rebin_is_charge_conserving_but_not_class_ready():
    result=run()
    assert result["CLASS_ready"] is False
    for item in result["rows"]:
        assert abs(item["absolute_number_error"])<1e-10
        assert item["min_fermion_occupation"]>=0
        assert item["max_fermion_occupation"]<=1
        assert item["numerical_fermion_entropy_increase"]>=-1e-10
