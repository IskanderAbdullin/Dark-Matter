"""Regression checks for near-degenerate decays and halo profile library."""
from pathlib import Path
import sys
import numpy as np
from scipy.integrate import quad
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'kinetics'))
from late_decay_bose_screen import transition
from late_decay_instant_equilibrium import thermal,equilibrium,C3,C5
from late_decay_line_gain import fixed_tau,crossing_tau
from late_decay_medium_screen import track
from two_state_halo_equilibrium import pl
from sparc_two_state_pilot import build,Library,K86

def test_recoil_probability():
    e=np.r_[0.,np.geomspace(1e-7,100,200)]
    q=np.array([.0001,.1,1.,5.])
    for b in [1e-12,1e-5,.1,1.,10.]:
        p=transition(q,b,e)
        assert np.all(p>=0) and np.max(abs(p.sum(axis=1)-1))<1e-13

def test_invariant_line_gain():
    e=np.array([0.,.1,.4,1.,2.,20.]);f=np.array([0.,7.,1e7,2e8,1e3])
    for r in [1e2,1e5,1e8]:
        a,b,_=fixed_tau(r,.001,1e-18,e,f,5e-8)
        assert abs(b/a-1)<1e-9
    r=1e5; a,_,_=fixed_tau(r,.001,1e-18,e,f,5e-8)
    b=crossing_tau(r,.001,1e-18,e,f,5e-8,order=24)
    assert abs(b/a-1)<1e-4

def test_bose_angular_integrals():
    for q,b,T,a in [(1.,.03,.7,0.),(.2,.3,2.,.02),(2.,1.,.1,.3),(.1,.1,3.,.001),(2.,1e-7,1.,.3)]:
        F,G=thermal(np.array([q]),b,T,a)
        f=lambda x: 1/np.expm1((q*q+b*b+2*q*b*x)/T+a)
        x=quad(f,-1,1,epsabs=1e-12)[0]/2
        y=quad(lambda x:(q*q+b*b+2*q*b*x)*f(x),-1,1,epsabs=1e-12)[0]/2
        assert np.isclose(F[0],x,rtol=3e-6,atol=1e-12)
        assert np.isclose(G[0],y,rtol=3e-6,atol=1e-12)

def test_ideal_BE_closure():
    for a0 in [.001,.1,1.,10.]:
        T0=2.3;n=C3*T0**1.5*pl(1.5,a0);K=C5*T0**2.5*pl(2.5,a0)
        T,a,Nc=equilibrium(n,K)
        assert Nc==0
        assert abs(C3*T**1.5*pl(1.5,a)/n-1)<2e-5
        assert abs(C5*T**2.5*pl(2.5,a)/K-1)<2e-5
    n=100.;K=1.;T,a,Nc=equilibrium(n,K)
    assert a==0 and Nc>0 and abs(C3*T**1.5*pl(1.5,a)+Nc-n)<1e-12

def test_density_cap_scaling():
    x=track(1e-5,.3);y=track(8e-5,.3)
    assert abs(y['r_unlock']/x['r_unlock']-2)<1e-6
    assert x['integrated_coherent_loss_fraction_no_backreaction']<1e-4
    assert x['locked_mixture_contact_sound_horizon_to_equality_Mpc']<=x['B_only_contact_sound_horizon_to_equality_Mpc']

def test_degenerate_two_state_halo_mapping():
    lib=Library(build(nr=300))
    a=.2;v=44.;m=86.;x=np.geomspace(.01,10,40)
    n1=lib.nu(a,0.,False);n2=lib.nu(a,0.,True)
    assert n2==2*n1
    r1=np.sqrt(K86/(v*n1))
    r2=np.sqrt(K86*(86/(m/2**.25))**4/(v*n2))
    assert abs(r1/r2-1)<1e-12
    assert np.all(lib.shape(a,0.,x)>0)


def test_collective_gain_limits():
    from late_decay_coherent_gain import dispersion
    d=np.array([1.,2.,3.]);w=np.array([1.,2.,1.])
    O=1e-3;s=dispersion(O,d,w)
    weak=np.pi*O*O*np.average(1/d,weights=w)/2
    assert abs(s/weak-1)<1e-5
    assert abs(dispersion(1e4,d,w)/1e4-1)<1e-7
