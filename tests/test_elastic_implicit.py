"""Stiff fixed-geometry conservative operator regression tests."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'kinetics'))
from elastic_implicit_trial import ImplicitElastic,evolve

def test_analytic_elastic_jacobian():
    m=ImplicitElastic(np.geomspace(1e-4,5,27),.02)
    r=m.validate_jac()
    assert r['jac_relative_error']<2e-5
    f=np.exp(-m.K)*1e10
    J=m.rates_jac(f,True)
    d=np.random.default_rng(42).normal(size=m.n)*f
    eps=2e-6
    fd=(m.rates_jac(f+eps*d,False)-m.rates_jac(f-eps*d,False))/(2*eps)
    assert np.linalg.norm(fd-J@d)/np.linalg.norm(fd)<1e-5
    v=m.vol
    assert np.max(abs(v@J)/np.maximum(v@abs(J),1e-300))<1e-12
    assert np.max(abs((v*m.K)@J)/np.maximum((v*m.K)@abs(J),1e-300))<1e-12

def test_implicit_short_evolution_conservation():
    e=np.geomspace(.02,10,25)
    q=np.sqrt(e[:-1]*e[1:]);f=np.exp(-q)*np.where(q<.5,50.,.1)
    r=evolve(e,f,.1,1.,.02,rtol=1e-6,wall=20.)
    assert r['success'],r['message']
    assert abs(r['number_relative_error'])<1e-10
    assert abs(r['energy_relative_error'])<1e-10
    assert r['maximum_negative_fraction']<1e-10
