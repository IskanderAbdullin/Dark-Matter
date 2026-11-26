"""Reproducible conditional massive-lepton kernel regression, not a full plasma."""
import math
import numpy as np
from scipy.special import expit
from kinetics.finite_mass_chie_kernel import kernel
from kinetics.microreversible_chie_lattice import invariant_kernel

def test_massless_limit_against_prior_kernel():
    for E,w,Ep in [(1.1,2.,1.7),(3.,4.,4.),(.7,1.2,1.4)]:
        wp=E+w-Ep
        a=float(kernel(.1,0.,.001,E,w,np.array([Ep]),np.array([wp]))[0])
        b=float(invariant_kernel(.1,.001,E,w,np.array([Ep]),np.array([wp]),nq=24)[0])
        assert math.isclose(a,b,rel_tol=2e-13)

def test_massive_reverse_and_pauli_equilibrium():
    M=.8;m=.00051099895;G=.001
    for E,w,Ep in [(1.2,.7,1.5),(2.4,1.3,1.8),(1.1,.02,1.105),(.801,.001,.8012)]:
        wp=E+w-Ep
        f=float(kernel(M,m,G,E,w,np.array([Ep]),np.array([wp]))[0])
        r=float(kernel(M,m,G,Ep,wp,np.array([E]),np.array([w]))[0])
        assert math.isclose(f,r,rel_tol=2e-12,abs_tol=1e-24)
        T=.4;mu_chi=.04;mu_e=-.08
        F=expit(-(E-mu_chi)/T)*expit(-(w-mu_e)/T)*expit((Ep-mu_chi)/T)*expit((wp-mu_e)/T)
        R=expit(-(Ep-mu_chi)/T)*expit(-(wp-mu_e)/T)*expit((E-mu_chi)/T)*expit((w-mu_e)/T)
        assert math.isclose(F,R,rel_tol=1e-13)

def test_forbidden_energy_transfer_and_threshold():
    M=.8;m=.00051099895;G=.001
    assert float(kernel(M,m,G,1.2,.7,np.array([1.89]),np.array([.01]))[0])==0.
    for eps in [1e-4,1e-6]:
        w=m+eps;E=1.2;Ep=E+.5*eps;wp=w-.5*eps
        assert float(kernel(M,m,G,E,w,np.array([Ep]),np.array([wp]))[0])>=0

def test_massless_continuity():
    E,w,Ep=1.1,2.,1.7;wp=E+w-Ep
    a=float(kernel(.1,0.,.001,E,w,np.array([Ep]),np.array([wp]))[0])
    b=float(kernel(.1,1e-7,.001,E,w,np.array([Ep]),np.array([wp]))[0])
    assert math.isclose(a,b,rel_tol=2e-12)
