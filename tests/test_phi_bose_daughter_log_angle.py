"""Independent exact angular controls for light real Bose decay daughters."""
import math
import pytest
from kinetics.phi_bose_daughter_log_angle import angular_nodes,analytic_angular_moments,weak_rate


@pytest.mark.parametrize('p',(.01,.1,1.,2.,8.,16.))
def test_log_angle_exact_rayleigh_jeans_moment_and_physical_positive_weights(p):
    r=analytic_angular_moments(p,.1,.02,86e-9,n=8)
    assert abs(r['RJ_1_over_E']/r['exact_RJ_1_over_E']-1)<1e-10
    assert abs(r['measure']/2-1)<2e-7
    c,w=angular_nodes(p,.1,.02,86e-9,n=8)
    assert min(c)>-1 and max(c)<1
    assert min(w)>0


def test_uniform_outgoing_cm_misses_bose_rayleigh_jeans_tail():
    r=analytic_angular_moments(16.,.1,.02,86e-9,n=8,scheme='uniform_cos')
    assert r['RJ_1_over_E']/r['exact_RJ_1_over_E']-1<-.5


def test_full_thermal_parent_rate_uniform_angular12_underestimates_ten_percent():
    kw=dict(pmax=2.,nr=48)
    old=weak_rate('thermal',nangle=12,scheme='uniform_cos',**kw)
    mapped12=weak_rate('thermal',nangle=12,scheme='log_ephi',**kw)
    mapped64=weak_rate('thermal',nangle=64,scheme='log_ephi',**kw)
    assert abs(mapped12['Rphi']/mapped64['Rphi']-1)<1e-7
    assert old['Rphi']/mapped64['Rphi']-1<-.09
    assert abs(old['H']/mapped64['H']-1)>.10
    assert mapped12['H']>0


def test_log_angular_does_not_cure_narrow_nonthermal_bose_shape_at_eight_nodes():
    kw=dict(pmax=16.,nr=48)
    a=weak_rate('sharp',nangle=8,scheme='log_ephi',**kw)
    b=weak_rate('sharp',nangle=96,scheme='log_ephi',**kw)
    c=weak_rate('sharp',nangle=48,scheme='log_ephi',**kw)
    assert abs(a['Rphi']/b['Rphi']-1)>.05
    assert abs(c['Rphi']/b['Rphi']-1)<1e-4
    assert b['H']>0


def test_parent_cutoff2_to16_remains_large_even_after_mapped_angle():
    a=weak_rate('sharp',pmax=2.,nr=48,nangle=96)
    b=weak_rate('sharp',pmax=16.,nr=48,nangle=96)
    assert abs(a['Rphi']/b['Rphi']-1)>.6
    assert a['Rphi']<0 and b['Rphi']<0
