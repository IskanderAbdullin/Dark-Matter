"""Nonlinear Zprime/finite ν thermal closure: charges, entropy and negative dt."""
import math
import pytest
from kinetics.vector_neutrino_finite_bath_time import integrate,coeff_thermal

def test_closed_finite_bath_number_energy_entropy_and_late_equilibrium():
    a=integrate(total_tau=2000,steps=200,n=64)
    b=integrate(total_tau=3000,steps=300,n=96)
    assert a["max_abs_N_relative_error"]<2e-12
    assert a["max_abs_E_relative_error"]<2e-12
    assert a["final"]["entropy_over_initial"]>1
    assert b["entropy_monotonic"]
    assert b["max_abs_N_relative_error"]<2e-12
    assert b["max_abs_E_relative_error"]<2e-12
    assert abs(b["delta_T_to_equilibrium"])<1e-5
    assert abs(b["delta_y_to_equilibrium"])<1e-5

def test_first_order_bath_splitting_convergence_and_negative_bias():
    a=integrate(total_tau=100,steps=20,n=64)
    b=integrate(total_tau=100,steps=40,n=64)
    c=integrate(total_tau=100,steps=80,n=64)
    ta,tb,tc=[v["final"]["T"] for v in (a,b,c)]
    assert abs(ta-tb)>1e-4
    assert 1.5<abs(ta-tb)/abs(tb-tc)<3

def test_oversized_split_step_fails_closed():
    with pytest.raises(ArithmeticError,match="step too large"):
        integrate(total_tau=3000,steps=20,n=64)

def test_nonzero_nu_number_chemical_potential_vs_independent_96_angular_gauss():
    import numpy as np
    from scipy.special import expit
    from numpy.polynomial.legendre import leggauss
    m=.03;T=1.;y=-1.1359
    p=np.array([.03,.3,1.,3.])
    a=coeff_thermal(m,p,T,y)
    c,w=leggauss(96)
    for pi,ai in zip(p,a):
        E=math.hypot(pi,m)
        E1=(E+pi*c)/2;E2=(E-pi*c)/2
        f1=expit(y-E1/T);f2=expit(y-E2/T)
        ref=float(w@(1-f1-f2)/2)
        assert math.isclose(ai,ref,rel_tol=2e-13)
