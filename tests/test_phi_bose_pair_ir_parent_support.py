"""Global physical massive real Bose daughter IR support for Dirac pair collisions."""
import math
import numpy as np
from kinetics.phi_bose_pair_ir_parent_support import (
    pair_rest_parent_threshold,global_pair_scalar_support_floor,
    pair_endpoint_scalar_shell,_generic_min_at_total_lab_energy_and_s
)
from kinetics.phi_bose_physical_cutcell_weak import M1,M2,MPHI

def test_exact_analytic_physical_pair_scalar_support_for_chi1_chi2_three_parent_cutoffs():
    for m,reference in ((M2, (.00009999751861499219,
                              .00002500003489762207,.000012500143037191316)),
                        (M1,(.0024984394508181726,
                             .0006249755908030623,.0003124969542185907))):
        for P,actual in zip((2.,8.,16.),reference):
            x=global_pair_scalar_support_floor(P,m)
            assert math.isclose(x,actual,rel_tol=3e-11)
            E=2*math.hypot(P,m)
            s=4*m*m
            check=_generic_min_at_total_lab_energy_and_s(E,s,MPHI)
            assert math.isclose(check,x,rel_tol=2e-8)
            assert math.isclose(pair_endpoint_scalar_shell(P,m),
                                x,rel_tol=3e-11)

def test_global_above_threshold_is_true_mass_not_false_parent_endpoint_min():
    for m,ref in ((M2,4651.162790654675),
                  (M1,116279.06976739886)):
        Pstar=pair_rest_parent_threshold(m)
        assert math.isclose(Pstar,ref,rel_tol=2e-14)
        assert math.isclose(global_pair_scalar_support_floor(Pstar,m),
                            MPHI,rel_tol=4e-13)
        for P in (Pstar*1.1,Pstar*2,1.5e5):
            assert math.isclose(global_pair_scalar_support_floor(P,m),
                                MPHI,rel_tol=3e-13)
            assert pair_endpoint_scalar_shell(P,m)>MPHI

def test_pair_chi2_opens_earlier_Bose_IR_than_restricted_decay_at_P2_to_P16():
    from kinetics.phi_bose_invariant_parent_support import global_support_floor
    for P in (2.,8.,16.):
        k2=global_pair_scalar_support_floor(P,M2)
        k1=global_pair_scalar_support_floor(P,M1)
        old=global_support_floor(MPHI,P)
        assert MPHI<k2<old<k1
        assert 11<old/k2<13
    assert pair_rest_parent_threshold(M2)<pair_rest_parent_threshold(M1)

def test_global_minimization_over_physical_E_and_s_matches_exact_bound():
    for m,P in ((M2,2.),(M2,16.),(M1,8.)):
        floor=global_pair_scalar_support_floor(P,m)
        # E total <=2 E(pmax), s >=4m². Physical subset of these
        # constraints is enough to prove a conservative global bound.
        Emax=2*math.hypot(P,m)
        for E in np.linspace(2*m,Emax,19):
            for s in np.linspace(4*m*m,E*E,23):
                x=_generic_min_at_total_lab_energy_and_s(E,s,MPHI)
                assert x>=floor-1e-10
    for m in (M2,M1):
        x=global_pair_scalar_support_floor(0.,m)
        assert math.isclose(x,m,rel_tol=3e-12)
