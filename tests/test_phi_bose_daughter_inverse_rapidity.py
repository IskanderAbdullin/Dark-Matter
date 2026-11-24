"""Independent massive χ2 two-sided parent rapidity and physical weak N/Q."""
import math
from kinetics.phi_bose_daughter_inverse_rapidity import (
    parent_bounds,daughter_coeff,exact_daughter_step,
    integrate_daughter_weak,integrate_daughter_adaptive,E2STAR,XI2,
)
from kinetics.phi_bose_physical_cutcell_weak import (
    M1,M2,KSTAR,analytic_weak,
)


def test_exact_finite_two_sided_parent_rapidity_support():
    pstar=M1*KSTAR/M2
    atrest=parent_bounds(M2,2.)
    assert math.isclose(atrest[0],pstar,rel_tol=1e-12)
    assert math.isclose(atrest[1],pstar,rel_tol=1e-12)
    for e in (.0201,.03,.052):
        lo,hi=parent_bounds(e,2.)
        assert 0<=lo<hi<2.
        if math.isclose(e,E2STAR,rel_tol=1e-12):
            assert lo<1e-9
        # Independently reconstruct Lorentz boosted on-shell energy shell.
        for p in (lo*.8,(lo+hi)/2,min(2.,hi+.1)):
            e1=math.hypot(p,M1)
            e_min=(e1*E2STAR-p*KSTAR)/M1
            e_max=(e1*E2STAR+p*KSTAR)/M1
            yes=e_min-1e-14<=e<=e_max+1e-14
            assert yes == (lo-1e-12<=p<=hi+1e-12)
    assert parent_bounds(20.,2.) is None


def test_conditional_exact_pauli_mode_step_and_chemical_balance():
    for e in (.03,.052,.1,.4):
        G,L,_,_=daughter_coeff(e,2.,48,T1=1.,Tphi=1.,
                                 mu1=-.15,muphi=-.2)
        f_eq=1/(1+math.exp(e-.05))
        assert G>0 and L>0
        assert math.isclose(G/(G+L),f_eq,rel_tol=2e-12)
        for dt in (0.,50.,1e6):
            f=exact_daughter_step(e,2.,dt,f_eq,48,
                                  T1=1.,Tphi=1.,mu1=-.15,muphi=-.2)
            assert 0<=f<=1
            assert math.isclose(f,f_eq,rel_tol=3e-12)
    for e in (.0201,.03,.052,.1,.4):
        for f0 in (0.,.3,1.):
            assert 0<=exact_daughter_step(e,2.,10000.,f0,32)<=1


def test_independent_adaptive_chi2_rapidity_first_all_three_energy_moments():
    for P in (2.,8.,16.):
        v=integrate_daughter_adaptive(P,36)
        a=analytic_weak(P)
        errs={k:abs(v[k]/a[k]-1)
              for k in ("Rphi","Qphi","Q1","Q2")}
        assert max(errs.values())<5e-9,(P,errs)
        assert abs(v["energy_defect"])<2e-19
        assert v["Rphi"]<0 and v["Qphi"]<0


def test_coarse_rapidity_gauss_misses_nontrivial_Bose_IR():
    for P in (2.,8.,16.):
        ref=analytic_weak(P)
        coarse=integrate_daughter_weak(P,8,8)
        fine=integrate_daughter_weak(P,64,64)
        e0=max(abs(coarse[k]/ref[k]-1)
               for k in ("Rphi","Qphi","Q1","Q2"))
        e1=max(abs(fine[k]/ref[k]-1)
               for k in ("Rphi","Qphi","Q1","Q2"))
        assert e1<e0/10
    r=integrate_daughter_weak(2.,64,64)
    ref=analytic_weak(2.)
    assert max(abs(r[k]/ref[k]-1) for k in
               ("Rphi","Qphi","Q1","Q2"))>1e-5
