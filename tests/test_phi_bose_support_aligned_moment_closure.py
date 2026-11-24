"""Support-aligned Bose DG-P1 fixed-fermion physical moment closure tests."""
import math
import numpy as np
from kinetics.phi_bose_moment_realisability_boundary import (
    boundary_crossing,dual_gap_at_fixed_right,
)
from kinetics.phi_bose_support_aligned_moment_closure import (
    fit_positive_etas,positive_eta_moments,target_moments,tangent_segment,
)
from kinetics.phi_bose_physical_cutcell_weak import MPHI,shell,cutcell_weak


def test_exact_support_edge_reverses_false_negative_first_eta_tangent():
    k0=shell(2.)[1]
    assert math.isclose(k0,.0011992509375794726,rel_tol=2e-13)
    for n in (24,32,40):
        v=cutcell_weak(2.,10,16,n)
        m,U=v['edges'][:2]
        R,Q=v['bycell'][0,:2]
        assert R<0 and Q<0
        _,_,old=tangent_segment(m,U,R,Q)
        C,G,new=tangent_segment(k0,U,R,Q)
        assert old[0]<0
        assert new[0]>0 and new[1]>0
        assert np.linalg.eigvalsh(G)[0]>0
        assert np.linalg.norm(G@new+C,ord=np.inf)<2e-25
    # Beneath kmin(P) this isolated 1->2 channel has zero on-shell events.
    assert MPHI<k0<cutcell_weak(2.,10,16,40)['edges'][1]


def test_finite_time_positive_two_moment_fits_100x_old_boundary():
    k0=shell(2.)[1]
    for n in (24,32,40):
        U=cutcell_weak(2.,10,16,n)['edges'][1]
        tstar=boundary_crossing(n)['boundary_time']
        for mul in (1.1,2.,10.,100.):
            v=fit_positive_etas(k0,U,tstar*mul)
            assert v['converged']
            assert v['fitted_eta_left']>k0 and v['fitted_eta_right']>U
            assert max(abs(x) for x in v['relative_moment_residual'])<5e-12
    assert dual_gap_at_fixed_right(
        32,1.1*boundary_crossing(32)['boundary_time'])['target_minus_max_left_M']>0


def test_support_aligned_N32_two_moment_fit_is_not_pointwise_exact_spectrum():
    k0=shell(2.)[1]
    U=cutcell_weak(2.,10,16,32)['edges'][1]
    t=2.*boundary_crossing(32)['boundary_time']
    fit=fit_positive_etas(k0,U,t)
    target=target_moments(k0,U,t)
    approx=positive_eta_moments(
        k0,U,fit['fitted_eta_left'],fit['fitted_eta_right'])
    assert np.allclose(approx,target,atol=1e-22,rtol=5e-12)
    assert math.isclose(
        fit['fitted_eta_left'],.0011994062360600238,rel_tol=2e-9)
    assert math.isclose(
        fit['fitted_eta_right'],.0021516295039263203,rel_tol=2e-9)
