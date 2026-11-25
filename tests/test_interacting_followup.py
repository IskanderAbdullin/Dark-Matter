"""Regression tests for the new conditional nonperturbative screens."""
from pathlib import Path
import sys
import numpy as np
from scipy.special import zeta
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'kinetics'))
from coherent_floquet import band,peak,floquet,oscillator
from radial_lattice_map import occupation_radial
from hartree_portal_pilot import HartreePortal
from halo_bose_rayleigh import quotient

def test_lame_quartic_band_and_monodromy():
    lo,hi=band(.5)
    assert np.isclose(lo,1.5) and np.isclose(hi,np.sqrt(3))
    p=peak(.5)
    assert abs(p['mu']-.0359794)<2e-7
    assert abs(p['det']-1)<1e-8
    assert floquet(.5,1.3)['mu']<1e-9
    assert floquet(.5,1.9)['mu']<1e-9

def test_harmonic_action():
    r=oscillator(.012,1.7,0.)
    assert np.isclose(r['action'],.012/1.7,rtol=1e-10)
    assert r['elliptic_parameter']==0

def test_radial_shell_number():
    n=32;L=24.;p0=.01
    k=2*np.pi*np.fft.fftfreq(n,d=L/n)
    k2=sum(a*a for a in np.meshgrid(k,k,k,indexing='ij'))
    e=np.array([0.,.01,.04,.08,.2]);f=np.array([0.,200.,800.,.1])
    fg,info=occupation_radial(k2,e,f,p0,L,6.)
    target=(200*(.04**3-.01**3)+800*(.06**3-.04**3))/(6*np.pi**2)
    assert np.isclose(fg.sum()*p0**3/L**3,target,rtol=1e-13)
    assert np.isclose(info['radial_represented_number'],target,rtol=1e-13)
    assert fg[0,0,0]==0 and np.all(fg>=0)

def test_hartree_stationary_energy_and_event_balance():
    m=HartreePortal(table=None,n=48,nparent=48,kappa=.00685,scalar_ratio=1e-12,
                   fermion_ratio=.9,qmin=1e-5,sigma=1e-7,parent_grid='balanced',
                   fermion_grid='balanced',parent_floor=1e-18)
    x=.2;p=1e8*np.exp(-(m.mean/.002)**2)
    y=np.r_[p,np.full(m.n+m.n1,.1),np.zeros(4)]
    assert m.check_local(x,y)['energy_gradient_relative_error']<1e-7
    d=m.dark(y[:m.dim])
    assert abs(m.energy_weights@d)/max(m.energy_weights@abs(d),1e-300)<1e-12

def test_halo_classical_limit_and_near_critical_convergence():
    r=quotient(1e-5,64)
    expected=(8/15)*np.sqrt(8/np.pi)
    assert abs(r['rayleigh_over_n_sigma_vT']/expected-1)<1e-5
    a=quotient(.99999999,64);b=quotient(.99999999,96)
    assert abs(a['rayleigh_over_n_sigma_vT']/b['rayleigh_over_n_sigma_vT']-1)<1e-5
    assert abs(b['n_integral']/(np.sqrt(np.pi)/2*zeta(1.5,1))-1)<2e-4
