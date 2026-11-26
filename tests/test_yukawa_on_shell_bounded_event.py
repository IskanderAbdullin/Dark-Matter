"""One elementary reaction conserves on-shell energy and real scalar Bose H."""
import math
import numpy as np
from scipy.special import expit
from kinetics.yukawa_on_shell_bounded_event import implicit_event,ledger,quantum_J
from kinetics.moving_yukawa_bose_pauli_weak import daughter_energies


def args():
 m1=.1;m2=.02;mphi=.001;p=.3;c=.2
 e2,eph=daughter_energies(p,c,m1,m2,mphi)
 return dict(E1=math.hypot(p,m1),E2=float(e2),Ephi=float(eph),
            w1=.1,w2=.15,wphi=.2,W=.05,dt=10.)


def test_on_shell_direct_weak_bose_event_preserves_two_charges_energy_entropy():
 a=args();start=(.6,.3,.1)
 old=ledger(*start,*(a[k] for k in ('E1','E2','Ephi','w1','w2','wphi')))
 f1,f2,b,x=implicit_event(*start,**a)
 new=ledger(f1,f2,b,*(a[k] for k in ('E1','E2','Ephi','w1','w2','wphi')))
 assert x>0
 assert 0<f1<1 and 0<f2<1 and b>0
 assert new['entropy']>old['entropy']
 for k in ('Nfermion','N1_plus_phi','energy'):
  assert abs(new[k]-old[k])<1e-14
 assert abs(x-a['dt']*a['W']*quantum_J(f1,f2,b))<1e-14


def test_inverse_can_consume_phi_and_no_negative_occupations():
 a=args();f1,f2,b,x=implicit_event(.05,.8,2.,**a)
 assert x<0 and b<2
 assert 0<=f1<=1 and 0<=f2<=1 and b>=0


def test_independent_chemical_mu_Bose_Fermi_stationary():
 a=args();T=.8;mu2=.02;muphi=-.01;mu1=mu2+muphi
 f1=expit((mu1-a['E1'])/T)
 f2=expit((mu2-a['E2'])/T)
 b=1/np.expm1((a['Ephi']-muphi)/T)
 f1p,f2p,bp,x=implicit_event(f1,f2,b,**a)
 assert max(abs(f1-f1p),abs(f2-f2p),abs(b-bp))<1e-13
 assert abs(x)<1e-13


def test_empty_bose_inverse_does_not_forbid_positive_bose_creation():
 a=args()
 f1,f2,b,x=implicit_event(.4,.1,0.,**a)
 assert x>0 and b>0 and 0<=f1<=1 and 0<=f2<=1


def test_250_deterministic_stiff_and_stimulated_events_never_clip_or_decrease_H():
 rng=np.random.default_rng(20261122)
 a=args()
 for _ in range(250):
  f1,f2=rng.uniform(.001,.999,2)
  b=10**rng.uniform(-5,4)
  w1,w2,wphi=10**rng.uniform(-3,2,3)
  W=10**rng.uniform(-3,2)
  dt=10**rng.uniform(-4,3)
  setting=dict(E1=a['E1'],E2=a['E2'],Ephi=a['Ephi'],
               w1=w1,w2=w2,wphi=wphi,W=W,dt=dt)
  old=ledger(f1,f2,b,a['E1'],a['E2'],a['Ephi'],w1,w2,wphi)
  F1,F2,B,x=implicit_event(f1,f2,b,**setting)
  now=ledger(F1,F2,B,a['E1'],a['E2'],a['Ephi'],w1,w2,wphi)
  assert -1e-13<=F1<=1+1e-13 and -1e-13<=F2<=1+1e-13 and B>=0
  assert now['entropy']>=old['entropy']-1e-12
  for k in ('Nfermion','N1_plus_phi','energy'):
   assert math.isclose(now[k],old[k],rel_tol=5e-14,abs_tol=5e-13)
