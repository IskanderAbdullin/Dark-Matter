"""Finite-volume moving Yukawa 1<->2+φ nonlinear Bose/Pauli, not cosmology."""
import math
import numpy as np
from scipy.special import expit
import pytest
from kinetics.yukawa_bose_pauli_bounded_time import (
    Config, density_cells, thermal_initial, rhs, ledger, evolve, euler_guard
)
from kinetics.moving_yukawa_bose_pauli_weak import weak_moments


def test_positive_physical_species_capacities_and_exact_radial_volume():
    c=Config(bins=32)
    for mass,g,e,w in ((c.m1,c.g1,c.e1,c.w1),
                       (c.m2,c.g2,c.e2,c.w2),
                       (c.mphi,c.gphi,c.ephi,c.wphi)):
        assert np.all(w>0)
        exact=g/(6*math.pi**2)*(e[-1]**2-mass**2)**1.5
        assert math.isclose(np.sum(w),exact,rel_tol=2e-14)


@pytest.mark.parametrize('bins',(16,32,64))
def test_three_species_strong_moments_identical_to_independent_moving_weak(bins):
    c=Config(bins=bins,nr=12,nangle=16)
    f1,f2,b=thermal_initial(c,T1=1.1)
    d1,d2,db,x=rhs(c,f1,f2,b)
    weak=weak_moments(lambda E:expit(-np.asarray(E)/1.1),
                      lambda E:expit(-np.asarray(E)),
                      lambda E:1/np.expm1(np.asarray(E)),
                      m1=c.m1,m2=c.m2,mphi=c.mphi,y=c.y,
                      pmax=c.pmax,nr=c.nr,nangle=c.nangle,g1=c.g1)
    for key1,key2 in [('Rphi','R_phi'),('Rchi1','R_chi1'),('Rchi2','R_chi2'),
                      ('Qchi1','Qchi1'),('Qchi2','Qchi2'),('Qphi','Qphi'),
                      ('entropy_event','entropy')]:
        assert math.isclose(x[key1],weak[key2],rel_tol=1e-10)
    assert abs(x['entropy_nodes']/x['entropy_event']-1)<2e-13
    assert x['entropy_event']>0
    assert abs(x['number_residual_chi'])<1e-24
    assert abs(x['number_residual_phi'])<1e-24
    assert abs(x['heat_residual'])<2e-24
    assert abs(np.dot(d1,c.w1)+np.dot(d2,c.w2))<1e-24
    assert abs(np.dot(d1,c.w1)+np.dot(db,c.wphi))<1e-24
    assert abs(np.dot(d1*c.e1,c.w1)+np.dot(d2*c.e2,c.w2)+np.dot(db*c.ephi,c.wphi))<1e-24


def test_full_nonzero_mu_fixed_point_and_bose_gain_loss():
    c=Config(bins=24,nr=8,nangle=12)
    f1,f2,b=thermal_initial(c,T1=.8,T2=.8,Tphi=.8,mu1=.03,mu2=.05,muphi=-.02)
    d1,d2,db,x=rhs(c,f1,f2,b)
    assert abs(x['Rphi'])<1e-23
    assert abs(x['entropy_event'])<1e-31
    assert max(np.max(abs(z)) for z in (d1,d2,db))<1e-19


def test_forward_bose_production_actual_time_preserves_all_three_species_bounds():
    c=Config(bins=32,nr=12,nangle=12)
    init=thermal_initial(c,T1=1.1)
    out=evolve(c,init,end_time=5e7,initial_dt=2e7,atol=1e-6)
    assert out['accepted']>10 and out['rejected']>0
    assert out['chi2_number_created']>0 and out['phi_number_created']>0
    assert math.isclose(out['chi2_number_created'],out['phi_number_created'],abs_tol=2e-14)
    assert max(abs(out[k]) for k in ('d_N1_N2','d_N1_Nphi','d_energy'))<2e-13
    assert out['min_step_entropy_gain']>=0 and out['entropy_gain']>0
    assert 0<out['min_fermion']<out['max_fermion']<1
    assert out['min_bose']>0
    assert out['max_local_full_vs_2half_spectral_error']<=1e-6


def test_undercool_mother_stimulated_inverse_destroys_boson_without_negative_entropy():
    c=Config(bins=24,nr=8,nangle=12)
    init=thermal_initial(c,T1=.9)
    _,_,_,x=rhs(c,*init)
    assert x['Rphi']<0 and x['entropy_event']>0
    out=evolve(c,init,end_time=1e7,initial_dt=1e7,atol=1e-5)
    assert out['phi_number_created']<0 and out['chi2_number_created']<0
    assert out['entropy_gain']>0
    assert max(abs(out[k]) for k in ('d_N1_N2','d_N1_Nphi','d_energy'))<2e-13


def test_zero_bose_endpoint_fail_closed_and_offshell_threshold():
    c=Config(bins=16,nr=6,nangle=8)
    state=list(thermal_initial(c))
    state[2]=np.zeros_like(state[2])
    with pytest.raises(ValueError):
        rhs(c,*state)
    with pytest.raises(ValueError):
        Config(m1=.1,m2=.02,mphi=.081)
