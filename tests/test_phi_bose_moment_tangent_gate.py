"""A negative φ IR source need not have feasible entropy-dual rest η tangent."""
import math
import numpy as np
from kinetics.phi_bose_moment_tangent_gate import tangent,two_hat_example
from kinetics.phi_bose_nonuniform_ir_grid import nodes,implicit_trial
from kinetics.phi_bose_moment_galerkin import Mode,physical_moments
from kinetics.phi_bose_implicit_moment_step import (
    source_with_jacobian,physically_weighted_yukawa_events,
)


def test_exact_two_overlapping_hat_counterexample_tangent_sign():
    x=two_hat_example(a=3,b=1,c=2,q=1)
    assert np.allclose(x['eta_rate'],[-.2,.6],rtol=0,atol=1e-14)
    assert np.allclose(x['eta_rate'],x['analytic_rate'],rtol=0,atol=1e-14)
    assert x['rest_eta_decreases_despite_zero_rest_hat_source']


def test_real_86_ev_tree_rate_has_feasibility_limit_despite_all_negative_phi_source():
    x=tangent(16,'power2')
    assert x['phi_rest_collisional_moment_source']==0
    assert x['all_yukawa_net_phi_is_destruction_no_real_BEC']
    assert x['positive_scalar_hat_sources']==0
    assert x['physical_event_entropy']>0
    assert -3.1e-9 < x['phi_rest_eta_time_derivative'] < -2.9e-9
    assert 28<x['linearized_max_allowed_time']<30
    assert x['linearized_boundary_limited_node']==0
    assert abs(x['total_energy_defect'])<1e-23
    assert abs(x['dark_fermion_number_defect'])<1e-23
    assert abs(x['parent_plus_phi_number_defect'])<1e-23
    assert implicit_trial(16,'power2',dt=20)['accepted']
    assert not implicit_trial(16,'power2',dt=30)['accepted']


def test_scalar_rest_mass_scales_linearized_critical_bound():
    a=tangent(16,'power2',mphi=86e-9)
    b=tangent(16,'power2',mphi=86e-7)
    ratio=b['linearized_max_allowed_time']/a['linearized_max_allowed_time']
    assert 95<ratio<105
    assert b['all_yukawa_net_phi_is_destruction_no_real_BEC']


def test_real_moment_susceptibility_recovers_phi_source_in_its_tangent():
    n=16;upper=math.hypot(2.,.1)+.2;m=86e-9
    modes=[Mode(mi,nodes(mi,upper,n),degeneracy=4,fermion=True)
           for mi in (.1,.02)] + [Mode(m,nodes(m,upper,n,'power2'))]
    eta=[modes[0].nodes/.9,modes[1].nodes,modes[2].nodes]
    events,_=physically_weighted_yukawa_events(pmax=2.,nr=12,nangle=12)
    C,_,_=source_with_jacobian(modes,eta,events)
    Cp=C[2*(n+1):]
    original=physical_moments(modes[2],eta[2],order=24)
    rate=np.linalg.solve(original['dM_deta'],Cp)
    candidate=eta[2]+rate
    assert min(candidate)>0
    actual=physical_moments(modes[2],candidate,order=24)
    assert np.max(abs((actual['M']-original['M'])-Cp)) < 1e-16
