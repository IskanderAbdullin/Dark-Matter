"""Focused physical-shell Bose DG-P1 tangent tests; not a finite-time solver."""
import math
import numpy as np
from kinetics.phi_bose_exact_cutcell_tangent import tangent,first_cell_susceptibility
from kinetics.phi_bose_physical_cutcell_weak import scalar_edges


def test_radial_angular_convergence_of_first_cell_physical_tangent():
    for n in (16,24,32,40):
        results=[tangent(n,2.,nr,na) for nr,na in
                 ((6,8),(10,16),(16,32),(32,48))]
        r0=results[0]['Rphi_first_GeV4']
        eta0=results[0]['deta_left_right'][0]
        for v in results[1:]:
            assert math.isclose(v['Rphi_first_GeV4'],r0,rel_tol=3e-11)
            assert math.isclose(v['deta_left_right'][0],eta0,rel_tol=3e-10)
            assert all(x>0 for x in v['eigenvalues_G'])
        assert r0<0
    assert tangent(16)['deta_left_right'][0]>0
    for n in (24,32,40):
        assert tangent(n)['deta_left_right'][0]<0
    assert math.isclose(tangent(32)['deta_left_right'][0],
                        -3.2640745039529283e-9,rel_tol=2e-10)
    assert math.isclose(
        tangent(32)['linear_tangent_positive_eta_limit_GeV_inverse'],
        26.34743780996741,rel_tol=2e-10)


def test_negative_inverse_fugacity_tangent_with_both_depleting_sources():
    x=tangent(32)
    assert x['source_hat_left_right'][0]<0
    assert x['source_hat_left_right'][1]<0
    assert x['deta_left_right'][0]<0 and x['deta_left_right'][1]>0
    G=np.array(x['G']);C=np.array(x['source_hat_left_right'])
    assert np.linalg.norm(G@np.array(x['deta_left_right'])+C,
                          ord=np.inf)<1e-25


def test_true_IR_first_cell_kinematically_forbidden_at_large_N():
    for n in (48,64):
        x=tangent(n)
        assert not x['first_cell_kinematically_accessible']
        assert x['Rphi_first_GeV4']==0
        assert x['Qphi_first_GeV5']==0
        assert all(v==0 for v in x['deta_left_right'])


def test_independent_high_temperature_massless_susceptibility_limit():
    upper=scalar_edges(2.,32)[1]
    G=first_cell_susceptibility(upper)
    high_T=upper/(2*math.pi**2)*np.array([[1/3,1/6],[1/6,1/3]])
    assert np.max(abs(G/high_T-1))<1.5e-4
    assert np.all(np.linalg.eigvalsh(G)>0)
