"""IR-aware Bose mesh gains and independent near-critical implicit failure."""
import math
from kinetics.phi_bose_nonuniform_ir_grid import nodes,moment_weak_rate,implicit_trial
from kinetics.phi_bose_ir_shape_gate import reference


def test_thermal_affine_bose_independent_of_nonuniform_grid():
    ref=reference('thermal')
    for scheme in ('uniform','power2','power3','two_zone'):
        r=moment_weak_rate('thermal',32,scheme)
        assert math.isclose(r['Rphi'],ref['R_phi'],rel_tol=1e-12)
        assert math.isclose(r['H'],ref['entropy'],rel_tol=1e-12)
        assert r['Q_residual']<1e-23


def test_power3_nonuniform_ir_shape_improves_over_uniform_at_3_grids():
    ref=reference('infrared_curved')
    for n in (32,64,128):
        uniform=moment_weak_rate('infrared_curved',n,'uniform')
        focused=moment_weak_rate('infrared_curved',n,'power3')
        eu=abs(uniform['Rphi']/ref['R_phi']-1)
        ef=abs(focused['Rphi']/ref['R_phi']-1)
        assert ef < .2*eu
        assert focused['H']>0 and focused['Rphi']<0


def test_nonuniform_two_zone_also_improves_broader_shape():
    ref=reference('broader')
    a=moment_weak_rate('broader',64,'uniform')
    b=moment_weak_rate('broader',64,'two_zone')
    assert abs(b['H']/ref['entropy']-1) < .1*abs(a['H']/ref['entropy']-1)


def test_near_critical_phi_implicit_inversion_fail_closed_not_fake_BEC():
    ok=implicit_trial(16,'uniform',dt=1e3)
    notok=implicit_trial(16,'power2',dt=1e3)
    assert ok['accepted'] and ok['H_gain']>0
    assert not notok['accepted']
    assert 'Newton cannot decrease residual' in notok['failure']
    assert notok['not_a_real_condensate']


def test_failed_ir_solver_has_no_positive_bose_microevent_production():
    # Distinguish discretization/manifold obstruction from real condensate.
    import numpy as np
    from kinetics.phi_bose_moment_galerkin import Mode
    from kinetics.phi_bose_implicit_moment_step import (
        physically_weighted_yukawa_events,source_with_jacobian,
    )
    import math
    n=16;upper=math.hypot(2.,.1)+.2
    modes=[Mode(m,nodes(m,upper,n),degeneracy=4,fermion=True)
           for m in (.1,.02)]
    modes.append(Mode(86e-9,nodes(86e-9,upper,n,'power2')))
    eta=[modes[0].nodes/.9,modes[1].nodes,modes[2].nodes]
    events,_=physically_weighted_yukawa_events(pmax=2.,nr=12,nangle=12)
    C,_,H=source_with_jacobian(modes,eta,events)
    phi=C[2*(n+1):]
    assert phi[0]==0
    assert all(x<=0 for x in phi)
    assert H>0 and sum(phi)<0


def test_nonuniform_critical_boundary_has_sharp_timestep_feasibility_limit():
    small=implicit_trial(16,'power2',dt=20.)
    too_large=implicit_trial(16,'power2',dt=30.)
    assert small['accepted'] and small['min_phi_eta']>0
    assert small['H_gain']>0 and max(abs(v) for v in small['N_E_defects'])<2e-12
    assert not too_large['accepted'] and too_large['not_a_real_condensate']
