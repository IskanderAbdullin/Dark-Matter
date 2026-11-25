"""Velocity-dependent resonance and late activation consistency tests."""
from pathlib import Path
import sys
import numpy as np
from scipy.integrate import quad
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'kinetics'))
from velocity_only_exposure import pair_intervals,pair_coefficient,halo_pdf
from resonant_mediator_inverse import width,averages,C,M,HBAR,SEC_GYR
from resonance_coupled_nr import Model
from direct_activation_consistency import sphere,direct

def test_pair_probability_and_maxwell():
    a,b,c=pair_intervals(np.array([1.]),np.array([1.]))
    x=np.linspace(0.,2.,2001)
    np.testing.assert_allclose(pair_coefficient(a,b,c,x[1:-1]),.5)
    assert abs(quad(lambda v:halo_pdf(v,1.),0,20)[0]-1)<1e-12

def test_resonance_width_roundtrip():
    p=width(90.)
    derivative=M*90e5/(2*C*C)
    area=np.pi*p['peak_unitarity_sigma_over_m']*p['width_eV']/(2*derivative)
    assert np.isclose(area,p['area'],rtol=1e-12)
    assert np.isclose(p['lifetime_Gyr']*SEC_GYR*p['width_eV'],HBAR)
    assert p['width_over_excess']<1e-20
    f,g=averages(np.array([.2]),.5,np.array([0.,1.]),np.array([5.]))
    np.testing.assert_allclose(f,25.)
    np.testing.assert_allclose(g,11.)

def test_coupled_resonance_invariants():
    m=Model(np.linspace(0,8,80),np.linspace(.001,12,52),1e6,1e6,1e-15,3)
    check=m.validate(0.,np.ones(m.na)*.1)
    assert check['balance_relative']<1e-10
    assert check['number_residual']<1e-12
    assert check['energy_residual']<1e-12
    assert check['jacobian_relative']<1e-5
    assert check['entropy_production']>=0

def test_direct_sphere_and_scaling():
    assert sphere(1e-7)==(1.,0.)
    assert abs(sphere(1e5)[0]*1e5-.75)<1e-10
    p=direct(.1,3.)
    q=direct(1.25,3.)
    assert np.isclose(q['vacuum_two_loop_threshold_ev4']/p['vacuum_two_loop_threshold_ev4'],np.sqrt(12.5))
    assert p['threshold_over_1percent_total']>1e5
