"""Numerical negative/positive controls for thermal Yukawa parent UV bound."""
import math
import pytest
from kinetics.phi_bose_chemical_uv_bound import forward_tail_upper, finite_tail_scan


def test_thermal_forward_bound_and_cutoff():
    rows=finite_tail_scan()
    assert len(rows)==6
    assert all(x['rate']>0 and x['tail_to_60']>0 for x in rows)
    assert all(x['strict_tail_bound']>=x['tail_to_60'] for x in rows)
    assert all(x['strict_tail_bound']>=x['tail_to_60']+
               x['finite_60_remainder_bound'] for x in rows)
    assert all(x['strict_tail_bound']>y['strict_tail_bound']
               for x,y in zip(rows,rows[1:]))
    assert rows[0]['relative_finite_60_tail']>.4
    assert rows[3]['relative_finite_60_tail']<2e-6
    assert rows[4]['upper_tail_fraction']<3e-7


def test_chemical_and_bose_domain_rejections():
    with pytest.raises(AssertionError):
        forward_tail_upper(8.,1.,-.21,.05,-.2)
    with pytest.raises(AssertionError):
        forward_tail_upper(8.,1.,-.09,.05,86e-9)
    with pytest.raises(AssertionError):
        forward_tail_upper(-1.,1.,-.09,.05,-.2)


def test_bound_is_not_finite_at_bose_critical_limit():
    # A finite-parent pmax gives a tighter kmin-based estimate, but the
    # universal b(mphi) bound intentionally diverges at mu_phi->mphi.
    vals=[]
    for eps in (1e-3,1e-4,1e-5):
        mup=86e-9-eps
        vals.append(forward_tail_upper(16.,1.,mup+.06,0.,mup))
    assert vals[0]<vals[1]<vals[2]
    assert vals[2]/vals[0]>50


def test_all_three_physical_energy_tail_certificates():
    from kinetics.phi_bose_chemical_uv_bound import forward_energy_tail_upper
    from kinetics.phi_bose_chemical_weak_reference import analytic
    full=analytic(60.,1.,1.,-.09,.05,-.2)
    for pmax in (2.,8.,16.):
        partial=analytic(pmax,1.,1.,-.09,.05,-.2)
        bound=forward_energy_tail_upper(pmax,1.,-.09,.05,-.2)
        for key in ('Qphi','Q1','Q2'):
            tail=full[key]-partial[key]
            assert abs(tail)<bound
        assert abs(sum(full[k]-partial[k] for k in ('Qphi','Q1','Q2')))<1e-23
