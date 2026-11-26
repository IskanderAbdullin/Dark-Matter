"""Implicit χ1χ2φ moment time, finite Bose/Pauli, entropy and energy tests."""
import math
import numpy as np
import pytest
from kinetics.phi_bose_implicit_moment_step import demonstration,implicit_step,source_with_jacobian,ledger
from kinetics.phi_bose_moment_galerkin import Mode,yukawa_moment_source


def system(bins=16,mphi=86e-9):
    modes=[Mode(m,np.linspace(m,4.,bins+1),fermion=i<2)
           for i,m in enumerate((.1,.02,mphi))]
    mus=(0.,.01,-.01)
    eta=[mo.nodes-mu for mo,mu in zip(modes,mus)]
    events=[(p,c,1e-4) for p in (0.,.05,.15,.7)
            for c in (-.8,-.2,.3,.9)]
    return modes,eta,events


def test_analytical_collision_jacobian_against_fd_columns():
    modes,eta,events=system(bins=8,mphi=.001)
    eta[0]=eta[0]+.12*np.sin(3*modes[0].nodes)
    C,J,H=source_with_jacobian(modes,eta,events)
    assert H>0
    for col in (0,10,24):
        lens=[len(m.nodes) for m in modes]
        x=np.concatenate(eta)
        ep=1e-5
        v=x.copy();v[col]+=ep
        vp=np.split(v,np.cumsum(lens)[:-1])
        v=x.copy();v[col]-=ep
        vm=np.split(v,np.cumsum(lens)[:-1])
        cp=source_with_jacobian(modes,vp,events)[0]
        cm=source_with_jacobian(modes,vm,events)[0]
        assert np.max(abs((cp-cm)/(2*ep)-J[:,col]))<2e-11


def test_physical_implicit_time_conserves_thermodynamic_charges_and_entropy():
    for dt in (.25,.5,1.,2.):
        run=demonstration(dt=dt)
        assert run['iterations']<20
        assert run['residual_abs']<5e-13
        assert abs(run['increments']['N1_plus_N2'])<5e-12
        assert abs(run['increments']['N1_plus_Nphi'])<5e-12
        assert abs(run['increments']['U_total'])<5e-12
        assert run['increments']['S_total']>0
        assert run['entropy_concavity_lower_bound']>0
        assert run['increments']['S_total']>=run['entropy_concavity_lower_bound']
        assert run['min_phi_eta']>0


def test_chemical_equilibrium_is_stationary_with_86_ev_phi():
    modes,eta,events=system(bins=16)
    r=implicit_step(modes,eta,events,10.,tol=1e-12)
    assert abs(r['increments']['N_phi'])<1e-12
    assert abs(r['increments']['S_total'])<1e-11
    assert r['residual_abs']<1e-12


def test_bose_critical_negative_inverse_fugacity_is_rejected_without_clipping():
    modes,eta,events=system(bins=8)
    eta[2][0]=-1e-4
    with pytest.raises(ValueError):
        implicit_step(modes,eta,events,1.)


def test_first_order_time_self_convergence_under_physical_bose_ir_moments():
    from kinetics.phi_bose_implicit_moment_step import time_grid_audit
    data=time_grid_audit(time_steps=(2,4,8,16,32,64))
    d=[x['Nphi_self_diff'] for x in data['self_differences']]
    assert all(x>0 for x in d)
    assert d[-2]/d[-1]>1.7
    assert d[-2]/d[-1]<2.3
    for row in data['rows']:
        assert row['min_step_entropy_gain']>=0
        assert row['min_final_phi_eta']>0
        assert max(abs(row['charge_N1N2_residual']),
                   abs(row['charge_N1Nphi_residual']),
                   abs(row['total_energy_residual']))<1e-11


def test_three_energy_grids_expose_nonconvergence_of_sparse_event_time_graph():
    from kinetics.phi_bose_implicit_moment_step import time_grid_audit
    vals=[time_grid_audit(bins=n,time_steps=(8,))['rows'][0]['delta_Nphi']
          for n in (8,16,32)]
    # The 16 isolated quadrature events are NOT dense under mode refinement;
    # a clean Gibbs number/H structure does not certify continuum evolution.
    assert abs(vals[1]/vals[2]-1)>.5
    assert all(v<0 for v in vals)


def test_physical_yukawa_parent_weight_matches_independent_weak_polar_integral():
    from kinetics.phi_bose_implicit_moment_step import physical_weak_rate_control
    for mphi in (86e-9,.001):
        z=physical_weak_rate_control(mphi=mphi,bins=16,nr=12,nangle=12)
        assert z['vacuum_gamma_rest']>0
        assert z['entropy_source']>0
        assert z['Rphi_from_moment_projection']<0
        assert abs(z['relative_discrepancy'])<5e-13


def test_against_preexisting_moving_yukawa_full_weak_reference():
    # This module exists unchanged in the actual Dark-Matter repository;
    # a small independent sandbox copy is used only for local testing.
    import numpy as np
    from scipy.special import expit
    from kinetics.moving_yukawa_bose_pauli_weak import weak_moments
    from kinetics.phi_bose_implicit_moment_step import physical_weak_rate_control
    f1=lambda E:expit(-np.asarray(E)/.9)
    f2=lambda E:expit(-np.asarray(E))
    b=lambda E:1/np.expm1(np.asarray(E))
    for mphi in (86e-9,.001):
        reference=weak_moments(f1,f2,b,m1=.1,m2=.02,mphi=mphi,y=.001,
                               pmax=2.,nr=12,nangle=12,g1=4)
        new=physical_weak_rate_control(mphi=mphi,bins=16,nr=12,nangle=12)
        assert math.isclose(new['Rphi_from_moment_projection'],
                            reference['R_phi'],rel_tol=2e-12)
        assert abs(reference['Qtotal'])<1e-24


def test_physically_weighted_time_uses_true_parent_width_but_is_not_mesh_converged():
    from kinetics.phi_bose_implicit_moment_step import physical_weight_time_audit
    x=physical_weight_time_audit(bins=(8,16,32),times=(1e7,))
    r=x['rows']
    assert len(r)==3 and x['vacuum_gamma_rest']>0
    assert all(row['Nphi_increment']<0 and row['S_increment']>0
               and row['min_phi_eta']>0 for row in r)
    assert all(max(abs(v) for v in row['charge_errors'])<5e-12
               and abs(row['energy_error'])<5e-12 for row in r)
    assert r[0]['Nphi_increment']/r[2]['Nphi_increment']>1.05
    # Different mode grids still alter the finite time result by many %.
    # This is a NEGATIVE physical-continuum gate, not a DM hypothesis test.
