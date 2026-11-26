"""Independent tests of physical Yukawa weak cut cells and IR aliasing.

All tests use fixed prescribed thermal states. They do not certify a nonlinear
strong solver, condensate, FRW or cosmological relic.
"""
import math
import numpy as np
from scipy.integrate import quad
from kinetics.phi_bose_physical_cutcell_weak import (
    analytic_weak, cutcell_weak, scalar_edges, shell, parent_splits, primitives,
)


def test_weak_N_and_all_three_Q_three_cutoffs_resolutions_reverse():
    settings=((.9,1.,0.,0.,0.), (1.,1.,-.09,.05,-.2),
              (1.,1.,-.21,.05,-.2))
    for s in settings:
        for pmax in (2.,8.,16.):
            ref=analytic_weak(pmax,*s)
            for nrad,nang,ncell in ((6,8,8),(10,16,16),(16,32,32)):
                v=cutcell_weak(pmax,nrad,nang,ncell,*s)
                for k in ('Rphi','Qphi','Q1','Q2'):
                    rel=abs(v[k]/ref[k]-1)
                    assert rel<(2e-7 if nrad==6 else 2e-9), (s,pmax,nrad,k,rel)
                assert v['min_event_weight']>0
                assert v['entropy_rate']>0
                assert abs(v['energy_defect'])<5e-24
                assert v['parent_segments']>=2
                assert v['capacities'].min()>0
                for k in range(4):
                    got=np.sum(v['bycell'][:,k])
                    assert math.isclose(got,v[('Rphi','Qphi','Q1','Q2')[k]],
                                        rel_tol=5e-14)


def test_inverse_and_nonzero_chemical_equilibrium():
    for pmax in (2.,8.,16.):
        fw=cutcell_weak(pmax,10,16,16,1.,1.,-.09,.05,-.2)
        rv=cutcell_weak(pmax,10,16,16,1.,1.,-.21,.05,-.2)
        eq=cutcell_weak(pmax,10,16,16,1.,1.,-.15,.05,-.2)
        assert fw['Rphi']>0 and fw['Qphi']>0
        assert rv['Rphi']<0 and rv['Qphi']<0
        assert fw['entropy_rate']>0 and rv['entropy_rate']>0
        assert max(abs(eq[k]) for k in ('Rphi','Qphi','Q1','Q2'))<2e-25
        assert abs(eq['entropy_rate'])<1e-35


def test_first_cell_exact_negative_control_and_mismatch_to_DOS_mean():
    v=cutcell_weak(2.,6,8,16)
    assert math.isclose(v['bycell'][0,0],-2.762149360447387e-12,rel_tol=2e-9)
    assert math.isclose(v['bycell'][0,1],-1.2277997619515968e-14,rel_tol=2e-9)
    assert math.isclose(v['Rphi'],-1.2954539554916677e-11,rel_tol=2e-9)
    assert len(parent_splits(2.,scalar_edges(2.,16)))==16
    # Constant-energy bin mean W/V cannot substitute for the on-shell
    # collision-weighted energy Qphi/Rphi without a real energy defect.
    effective=v['bycell'][0,1]/v['bycell'][0,0]
    dosmean=v['capacities'][0,1]/v['capacities'][0,0]
    assert dosmean/effective>1.45


def test_true_hard_parent_support_and_independent_daughter_integrals():
    v=cutcell_weak(2.,10,16,64)
    assert v['bycell'][0,0]==0 and v['bycell'][0,1]==0
    for p in (.1,1.,2.,8.,16.):
        E,lo,hi,_=shell(p)
        n,q=primitives(E,lo,hi,1.,.05,-.2)
        fd=lambda k:1/(1+math.exp(E-k-.05))
        bose=lambda k:1/math.expm1(k+.2)
        ni=quad(lambda k:fd(k)*bose(k),lo,hi,epsabs=1e-17,epsrel=2e-11)[0]
        qi=quad(lambda k:k*fd(k)*bose(k),lo,hi,epsabs=1e-17,epsrel=2e-11)[0]
        assert math.isclose(n,ni,rel_tol=3e-10)
        assert math.isclose(q,qi,rel_tol=3e-10)


def test_physical_parent_partition_eliminates_false_small_first_cell():
    naive=cutcell_weak(2.,48,16,16,split_parent=False)
    split=cutcell_weak(2.,6,8,16,split_parent=True)
    target=-2.762149360447387e-12
    assert abs(naive['bycell'][0,0]/target-1)>1e-5
    assert abs(split['bycell'][0,0]/target-1)<2e-9
