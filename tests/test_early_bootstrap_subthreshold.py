"""Numerical checks for the early occupancy and subthreshold screens."""
from pathlib import Path
import sys
import numpy as np
from scipy.integrate import quad
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'kinetics'))
from early_bose_bootstrap import rearrange
from early_hartree_coherence import equilibrium
from subthreshold_transport_tradeoff import angle_integral,rayleigh_grid,mb_integral,bose_source
from halo_bose_rayleigh import quotient
from soft_mixing_activation import scales,frac_B
from two_state_halo_equilibrium import source

def test_phase_space_rearrangement():
    for p in [.001,.1,1.,100.,1e5]:
        F=.01;m=86.;n=F*p**3/(6*np.pi**2)
        I,V=rearrange(n,F,m)
        refI=quad(lambda q:F*q/np.hypot(m,q),0,p,epsabs=1e-20)[0]
        refV=quad(lambda q:F*q*q/(2*np.pi**2*np.hypot(m,q)),0,p,epsabs=1e-20)[0]
        assert np.isclose(I,refI,rtol=1e-10,atol=0.)
        assert V>=refV*(1-1e-12)
        assert V<=n/m*(1+1e-12)

def test_hartree_fixed_point():
    p=scales(.3,z_on=8.)
    E=np.array([1000.,2000.]);w=np.array([.4,.6])
    D=p['Delta'];off=p['delta']*1e-8
    ds,u=equilibrium(E,w,1e14,D,off,p['lambda_B'])
    assert ds>=D and u<=frac_B(1.,off,D)
    assert np.isclose(u,frac_B(1.,off,ds),rtol=1e-12)

def test_subthreshold_angular_integral():
    for a in [1e-6,.001,.1,1.,1000.]:
        for l,h in [(0.,1.9),(.1,.1001),(1.999,1.99999)]:
            got=angle_integral(2.,l,h,a)
            ref=quad(lambda P:(a/(a+4-P*P))**2,l,h,epsabs=1e-22,epsrel=1e-10)[0]
            assert np.isclose(got,ref,rtol=3e-9,atol=0.)

def test_transport_contact_limit():
    Q,contact=rayleigh_grid(np.array([1.,1e7]),order=48)
    r=quotient(np.exp(-1e-6),48)['rayleigh_over_n_sigma_vT']
    assert np.isclose(contact,r,rtol=1e-9)
    assert np.isclose(Q[1],r,rtol=2e-6)
    assert Q[0]<Q[1]

def test_mb_threshold_and_bose_contact():
    for d in [.02,.4,5.,15.]:
        a=5.
        ref=quad(lambda y:2*np.sqrt(2/np.pi)*np.sqrt(y*(y+d))*np.exp(-y-d)*(a/(a+2*(y+d)))**2,0,80,epsabs=1e-20,epsrel=1e-9)[0]
        assert np.isclose(mb_integral(d,a,128),ref,rtol=.002,atol=0.)
    q=bose_source(.04,.3169873,1e7,48,48)
    ref=source(.04,.3169873,48,48,alpha=1e-6)['heavy_initial_number_gain_over_n_sigma_vT']
    assert np.isclose(q,ref,rtol=2e-6)
