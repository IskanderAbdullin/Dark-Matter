"""Physical upper-edge true Bose eta-dot versus false affine DG-P1 endpoint sign."""
import math
from kinetics.phi_bose_upper_support_curvature import (
    audit,curvature_fit,modal_eta_tangent,physical_weak_segment,
)
from kinetics.phi_bose_physical_cutcell_weak import cutcell_weak,shell


def test_true_modal_eta_dot_positive_and_vanishes_at_upper_support():
    K=shell(2.)[2]
    for d in (.1,.05,.025,.0125,.00625,.003125,.0015625):
        assert modal_eta_tangent(K-d)>0
    assert modal_eta_tangent(K)<1e-20
    A,B=curvature_fit()
    assert 2.38e-11<A<2.39e-11 and 7.88e-12<B<7.90e-12


def test_scalar_first_upper_supported_cell_matches_parent_first_physical_N_Q():
    V=cutcell_weak(2.,10,16,32)
    j=29
    l=V['edges'][j];K=shell(2.)[2]
    R,Q=physical_weak_segment(l,K,2.,48,48)
    assert math.isclose(R,V['bycell'][j,0],rel_tol=1e-11)
    assert math.isclose(Q,V['bycell'][j,1],rel_tol=1e-11)


def test_negative_right_eta_tangent_survives_alignment_and_refinement():
    d=audit()
    assert d['last_active_bin']==29
    rows=d['rows']
    for r in rows:
        assert r['eta_left_dot_GeV']>0
        assert r['eta_right_dot_GeV']<0
        assert r['Rphi_GeV4']<0
    for x,y in zip(rows[:-1],rows[1:]):
        assert abs(y['eta_right_dot_GeV'])<abs(x['eta_right_dot_GeV'])/3.8
    assert abs(rows[-1]['ratio_to_asymptotic']-1)<.002
    assert math.isclose(rows[5]['eta_right_dot_GeV'],
                        -1.6225765066796667e-17,rel_tol=1e-5)
