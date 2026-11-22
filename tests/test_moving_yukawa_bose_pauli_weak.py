"""Fully boosted Yukawa quantum event is not full φ cosmology."""
import math
import numpy as np
from scipy.special import expit
import pytest
from kinetics.moving_yukawa_bose_pauli_weak import vacuum_width,daughter_energies,event,weak_moments


def test_scalar_yukawa_threshold_and_exact_massless_width():
 assert vacuum_width(.1,.02,.081,.001)==0
 assert math.isclose(vacuum_width(.1,0,0,.001),1e-6*.1/(16*math.pi),rel_tol=1e-13)
 for p in (0.,.01,1.):
  e2,eph=daughter_energies(p,np.linspace(-1,1,17),.1,.02,.001)
  assert max(abs(e2+eph-math.hypot(.1,p)))<1e-13
  assert min(e2)>=.02 and min(eph)>=.001


def test_quantum_chemical_detailed_balance_nonneutral_mu_and_bose():
 T=.8;mu1=.03;mu2=.05;muphi=-.02
 f1=lambda E:expit((mu1-np.asarray(E))/T)
 f2=lambda E:expit((mu2-np.asarray(E))/T)
 bp=lambda E:1/np.expm1((np.asarray(E)-muphi)/T)
 for p in (0.,.2,2.):
  x=event(p,np.linspace(-.9,.9,11),.1,.02,.001,f1,f2,bp)
  assert max(abs(x['flux']))<1e-13
  assert max(abs(x['entropy']))<1e-25


def test_out_of_equilibrium_energy_and_stoichiometric_weak_moments():
 f1=lambda E:expit(-np.asarray(E)/1.1)
 f2=lambda E:expit(-np.asarray(E))
 bp=lambda E:1/np.expm1(np.asarray(E))
 x=weak_moments(f1,f2,bp,nr=16,nangle=16)
 assert x['gamma_rest']>0 and x['R_phi']>0
 assert x['R_chi1']==-x['R_chi2']==-x['R_phi']
 assert abs(x['Qtotal'])<1e-16
 assert x['entropy']>0
 assert x['max_event_energy_defect']<1e-13


def test_vacuum_bose_zero_inverse_and_bose_stimulation_gain():
 f1=lambda E:expit(-np.asarray(E)/.9)
 f2=lambda E:expit(-np.asarray(E))
 zero=lambda E:np.zeros_like(np.asarray(E))
 finite=lambda E:1/np.expm1(np.asarray(E))
 x=event(.5,np.linspace(-.8,.8,13),.1,.02,.001,f1,f2,zero)
 assert min(x['flux'])>0 and max(x['reverse'])==0
 assert np.all(np.isinf(x['entropy']))
 b=event(.5,np.linspace(-.8,.8,13),.1,.02,.001,f1,f2,finite)
 assert np.all(b['forward']>x['forward'])
