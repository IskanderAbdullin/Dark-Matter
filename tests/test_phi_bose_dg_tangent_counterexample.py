"""Nonoverlapping physical Bose cells still have within-cell forbidden tangent."""
import numpy as np
from kinetics.phi_bose_dg_tangent_counterexample import local_spd_counterexample,physical_case
from kinetics.phi_bose_cell_entropy import fixture,source,mode_moments
from kinetics.phi_bose_daughter_log_angle import physical_events


def test_exact_spd_two_local_hat_negative_rest_tangent():
    x=local_spd_counterexample(a=3,b=1,c=2,u=.05,v=1)
    assert x['C']==[-.05,-1.]
    assert np.allclose(x['eta_rate'],[-.18,.59],rtol=0,atol=1e-14)
    assert abs(x['analytical_first_eta_rate']+.18)<1e-14
    assert x['rest_eta_decreases_despite_both_local_sources_nonpositive']


def test_physical_mapped_n32_finer_mesh_fails_while_n16_passes():
    a=physical_case(16,48,run_time_tests=True)
    b=physical_case(32,48,run_time_tests=True)
    assert a['first_eta_derivative']>0
    assert b['first_eta_derivative']<0
    assert b['all_scalar_source_negative'] and b['H']>0 and b['source_Rphi']<0
    assert all(x<0 for x in b['first_cell_C'])
    assert abs(b['first_eta_derivative']-b['analytic_eta_derivative'])<1e-22
    assert 24<b['boundary_linear_dt']<25
    assert a['implicit_trials']['30.0']['accepted']
    assert b['implicit_trials']['20.0']['accepted']
    assert not b['implicit_trials']['30.0']['accepted']
    assert not b['implicit_trials']['1000.0']['accepted']


def test_angular_refinement_validates_mapped_rate_but_not_feasibility():
    for na in (12,24,48,96):
        x=physical_case(32,na,nrad=24)
        assert x['all_scalar_source_negative']
        assert x['first_eta_derivative']<0
        assert 24<x['boundary_linear_dt']<30
        assert x['H']>0
        assert abs(x['source_Rphi']/(-1.295453955424128e-11)-1)<2e-6


def test_uniform_mapped_angular_discrepancy_exposes_false_first_cell_safety():
    fake=physical_case(32,12,nrad=24,scheme='uniform_cos',run_time_tests=True)
    real=physical_case(32,48,nrad=24,scheme='log_ephi')
    assert fake['first_cell_C']==[0.,0.]
    assert fake['implicit_trials']['30.0']['accepted']
    assert real['first_eta_derivative']<0
    assert abs(fake['source_Rphi']/real['source_Rphi']-1)>.09


def test_tangent_mass_matrix_reproduces_negative_local_source_without_clipping():
    modes,eta,_=fixture(n=32)
    events,_=physical_events(pmax=2.,nr=24,nangle=48)
    C,_,H=source(modes,eta,events)
    n=modes[2].size
    cp=C[-n:]
    old=mode_moments(modes[2],eta[2],order=24)
    tangent=np.linalg.solve(old['J'],cp)
    assert H>0 and tangent[0]<0
    eps=.01
    after=mode_moments(modes[2],eta[2]+eps*tangent,order=24)
    assert np.max(abs(after['M']-old['M']-eps*cp))<1e-16


def test_nonmonotone_mesh_width_limits_negative_rest_tangent():
    grid=[(n,physical_case(n,48)['first_eta_derivative'])
          for n in (16,20,24,32,40,48,64)]
    assert grid[0][1]>0 and grid[1][1]>0
    assert all(x[1]<0 for x in grid[2:5])
    assert grid[5][1]==0 and grid[6][1]==0
