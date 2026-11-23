"""Independent regression tests for massive-boson inverse-decay chemical source."""
import math
from scipy.integrate import quad
from scipy.special import expit
from kinetics.phi_bose_chemical_weak_reference import analytic, gauss, kinematics, moments


def test_exact_primitives_against_direct_daughter_quadrature():
    for p in (.1,1.,2.,8.,16.):
        e,lo,hi,_=kinematics(p)
        i0,i1=moments(e,lo,hi,1.,.05,-.2)
        direct=lambda k: expit(-((e-k)-.05))/(math.expm1(k+.2))
        j0=quad(direct,lo,hi,epsabs=1e-17,epsrel=1e-11)[0]
        j1=quad(lambda k:k*direct(k),lo,hi,epsabs=1e-17,epsrel=1e-11)[0]
        assert math.isclose(i0,j0,rel_tol=3e-10)
        assert math.isclose(i1,j1,rel_tol=3e-10)


def test_chemical_equilibrium_reverse_sign_and_positive_entropy():
    eq=gauss(2.,24,12,1.,1.,-.15,.05,-.2)
    fwd=gauss(2.,24,12,1.,1.,-.09,.05,-.2)
    rev=gauss(2.,24,12,1.,1.,-.21,.05,-.2)
    assert abs(eq['Rphi'])<2e-25
    assert fwd['Rphi']>0 and rev['Rphi']<0
    for x in (eq,fwd,rev):
        assert abs(x['energy_residual'])<1e-23
        assert x['entropy_rate']>=-1e-24


def test_three_parent_cutoffs_and_three_radial_angular_resolutions():
    for pmax in (2.,8.,16.):
        ref=analytic(pmax,1.,1.,-.09,.05,-.2)
        errs=[]
        for n,a in ((24,12),(48,24),(72,48)):
            got=gauss(pmax,n,a,1.,1.,-.09,.05,-.2)
            errs.append(max(abs(got[k]/ref[k]-1)
                            for k in ('Rphi','Qphi','Q1','Q2')))
            assert abs(got['energy_residual'])<1e-23
            assert got['entropy_rate']>0
        assert errs[-1]<5e-9
        assert errs[-1]<errs[0]
