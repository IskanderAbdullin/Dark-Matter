"""Independent checks for softly mixed late activation and derivative portal."""
from pathlib import Path
import sys
import numpy as np
import mpmath as mp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'kinetics'))
from derivative_portal_consistency import widths,traces,harmonic_yield,cosine_yield
from derivative_forward_screen import check_trace,angle_amplitude
from soft_mixing_activation import scales,frac_B,bath_integrals
from soft_mixing_energy_check import bath_energy_table
from soft_mixing_stability_screen import thresholds
from two_state_halo_equilibrium import pl,box_equilibrium,source
from late_activation_initial_conditions import action

def test_derivative_vacuum_ward_traces():
    for m in [0.,.01,.05,.099]:
        np.testing.assert_allclose(widths(1.,.9,m,2.),traces(1.,.9,m,2.),rtol=1e-11)
    assert check_trace()<1e-12

def test_forward_soft_limit():
    k=1e-7;p=.4;mi=1.;mj=.9;F=2.
    D=mi**2-mj**2
    expect=-8*k*k*(mi*mi+4*p*p/3)/(D*F**2)
    assert np.isclose(angle_amplitude(k,p,mi,mj,F),expect,rtol=1e-10,atol=0.)

def test_soft_mixing_eigenvalue_and_threshold():
    p=scales()
    assert np.isclose(frac_B(1.,p['delta'],p['Delta']),p['u'])
    assert np.isclose(p['lambda_B']*p['u']**2,p['lambda_late'])
    mat=np.array([[0.,p['delta']],[p['delta'],p['Delta']]])
    ev,vec=np.linalg.eigh(mat)
    assert np.isclose(vec[1,0]**2,p['u'])
    assert np.isclose((ev[1]-ev[0])/(2*86),p['mass_gap_eV'],rtol=1e-8)
    t=thresholds(p)
    assert np.isclose(t['cutoff_over_m'],p['overlap_ratio']**.25)

def test_collisionless_eta_energy_derivative():
    v=bath_integrals();f=bath_energy_table()
    assert np.isclose(v(0.),1/12)
    for j in [.01,.1,1.,10.,100.]:
        d,vr=f(j);eps=j*1e-4
        der=(f(j+eps)[0]-f(j-eps)[0])/(2*eps)
        assert np.isclose(der,j*vr,rtol=2e-5)
        assert np.isclose(v(j),vr,rtol=2e-5)

def test_polylog_and_two_state_equilibrium():
    for s in [1.5,2.5]:
        for a in [0.,.01,.1,.49,.5,1.,10.]:
            assert np.isclose(pl(s,a),float(mp.polylog(s,mp.exp(-a))),rtol=2e-13)
    x=box_equilibrium(.02)
    assert abs(x['number_relative_residual'])<1e-9
    assert abs(x['energy_relative_residual'])<1e-9
    assert x['condensate_number_fraction']==0.
    assert .45<x['heavy_number_fraction']<.5

def test_halo_conversion_quadrature():
    a=source(.02,.3169873,48,48)
    b=source(.02,.3169873,96,96)
    assert abs(a['heavy_initial_number_gain_over_n_sigma_vT']/b['heavy_initial_number_gain_over_n_sigma_vT']-1)<.01

def test_spectator_harmonic_action():
    assert np.isclose(action(.4,1.,0.),.4,rtol=1e-11)
    assert np.isclose(action(.5,.25,0.),1.,rtol=1e-11)
