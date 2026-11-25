"""Checks for subthreshold scattering, three-body kinematics, and trigger energy."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'kinetics'))
from subthreshold_cluster_screen import cross_sections,SCALE
from subthreshold_recombination import final_state,norm2,amplitude,zero_rate,radial_bose
from direct_subthreshold_activation import sunset_threshold
from direct_subthreshold_background import energy
from resonant_mediator_inverse import M
from soft_mixing_activation import frac_B
from subthreshold_roll_scan import params

def test_subthreshold_cross_sections():
    u=1e-5;D=.006;delta=D*np.sqrt(u*(1-u))/(1-2*u);a=5.;sig=.4
    v=1000.
    r=cross_sections(sig,u,D,delta,a,v)
    # Detailed balance for opposite channels at equal invariant energy:
    # high v and tiny splitting imply the final weight ratios up to phase space.
    assert r['00->11']/r['00->00']>1e9
    low=cross_sections(sig,u,D,delta,a,1e-2)
    assert abs(low['00->00']/sig-1)<1e-6
    assert low['00->11']==0.
    assert low['00->01']==0.
    assert r['11->11']>r['00->11']
    assert all(np.isfinite(x) and x>=0 for x in r.values())

def test_three_body_kinematics_and_tree_symmetry():
    rng=np.random.default_rng(223);p=rng.normal(size=(3,8,3))*.01
    a,b,c=p;ms=np.sqrt(4*M*M-.001)
    vs,S,k,root,ok,E,po,ps=final_state(a,b,c,M,M,ms,np.linspace(-1.,1.,8),np.ones(8))
    assert np.all(ok)
    assert np.max(abs(sum(vs)+S))<1e-11
    for v in vs: assert np.max(abs(norm2(v)-M*M))<1e-8
    assert np.max(abs(norm2(S)-ms*ms))<1e-8
    g=.01;lam=3*g*g/ms**2
    A=amplitude(vs,S,1.,g,lam,M,M+1,ms)
    B=amplitude([vs[2],vs[0],vs[1],vs[3]],S,1.,g,lam,M,M+1,ms)
    np.testing.assert_allclose(A,B,rtol=1e-13)
    Z=zero_rate(1.,g,lam,M,M+1,ms)
    assert Z['allowed'] and Z['coefficient']>0
    Z2=zero_rate(1.,2*g,4*lam,M,M+1,ms)
    np.testing.assert_allclose(Z2['coefficient']/Z['coefficient'],64.,rtol=1e-12)

def test_sunset_at_equal_threshold():
    v=sunset_threshold(1.,M)
    np.testing.assert_allclose(v['I_eV2']/M**2,8*np.log(4)-15,rtol=1e-13)
    ref=M**2*(.25-(8*np.log(4)-15))/4/(16*np.pi**2)**2
    np.testing.assert_allclose(v['total_eV4'],ref,rtol=1e-13)
    np.testing.assert_allclose(sunset_threshold(2.,M)['total_eV4'],2*ref,rtol=1e-13)

def test_eta_radiation_continuity():
    p=params(10,.1,.7,1.,10)
    t=np.linspace(-np.log(101.),0.,4001);r=np.exp(-t)
    # A stationary x=0 field at the symmetric point has constant vacuum
    # energy and massless collisionless eta radiation.
    rr=dict(t=t,r=r,x=np.zeros_like(t),dx=np.zeros_like(t))
    v=energy(rr,p,'vev-curvature')
    assert v['max_abs_comoving_continuity_residual']<3e-6
    assert v['today_eta_fraction']>0

def test_bose_radial_normalization():
    _,_,meta=radial_bose(1e-6)
    assert abs(meta['number_rel_error'])<1e-5
    assert abs(meta['mean_q2_rel_error'])<1e-5
