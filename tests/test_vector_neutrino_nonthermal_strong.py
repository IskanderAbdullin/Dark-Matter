"""Independent nonthermal strong ν and weak Z' collision moment closure."""
import math
import numpy as np
import pytest
from scipy.special import expit
from kinetics.vector_neutrino_nonthermal_strong import (
    vector_weak,strong_nu_moments,strong_nu_mode,gl,GV,GNU)

@pytest.mark.parametrize("sign,r",[(1,.75),(-1,1.25),(1,1.),(0,.75)])
def test_independent_nonthermal_number_energy_entropy(sign,r):
    m=.3
    f=lambda k:expit(-k+.15*sign*np.exp(-k))
    b=lambda p:r/np.expm1(np.hypot(p,m))
    w=vector_weak(m,f,b,np_=64,nangle=64)
    s=strong_nu_moments(m,f,b,np_=64,nk=64)
    assert w["entropy_weak"]>0
    assert math.isclose(s["dNnu"],w["dNnu_expected"],rel_tol=1e-9)
    assert math.isclose(s["dEnu"],w["dEnu_expected"],rel_tol=1e-9)
    assert math.isclose(s["dSnu"]+w["entropy_vector"],w["entropy_weak"],rel_tol=1e-9)

def test_exact_inverse_shell_support():
    m=.3
    for p in (.01,.1,.5,2,10):
        E=math.hypot(p,m)
        for k in ((E-p)/2,(E+p)/2):
            assert math.isclose(abs(k-m*m/(4*k)),p,rel_tol=5e-12)
    assert (GV,GNU)==(3,4)

def test_unsplit_gauss_falsely_breaks_weak_strong_equality():
    m=.3
    f=lambda k:expit(-k+.15*np.exp(-k))
    b=lambda p:.75/np.expm1(np.hypot(p,m))
    ref=vector_weak(m,f,b,np_=64,nangle=64)["dNnu_expected"]
    pmax=12; Emax=math.hypot(pmax,m)
    ks,wk=gl((Emax-pmax)/2,(Emax+pmax)/2,32)
    unsplit=sum(GNU/(2*math.pi**2)*wi*k*k*
        strong_nu_mode(float(k),m,f,b,np_=64) for k,wi in zip(ks,wk))
    split=strong_nu_moments(m,f,b,np_=64,nk=32)["dNnu"]
    assert abs(unsplit/ref-1)>1e-4
    assert abs(split/ref-1)<1e-6
