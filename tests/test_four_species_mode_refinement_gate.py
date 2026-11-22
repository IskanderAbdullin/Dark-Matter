"""Energy mode, not 5D phase-space, interpolation consistency gate."""
import numpy as np
import pytest
from kinetics.four_species_mode_refinement_gate import (
    linear_interpolation_bound,evaluate_rate_refinement,
)


def test_analytic_Fermi_logit_interpolation_remainder_for_nonlinear_modes():
    for lower,amp,k,T in ((.1,.3,4.,1.1),(.02,.2,3.,.9),(0.,.25,4.,1.)):
        for n in (64,256,1024):
            gr=np.linspace(lower,40.1,n+1)
            h=gr[1]-gr[0]
            x=np.linspace(lower,40.1,10001)
            f=lambda E:-np.asarray(E)/T+amp*np.sin(k*np.asarray(E))
            err=np.max(abs(f(x)-np.interp(x,gr,f(gr))))
            assert err<=linear_interpolation_bound(h,amp,k)*(1+1e-8)


def test_nonthermal_energy_entropy_order_emerges_only_on_refined_mesh():
    x=evaluate_rate_refinement((256,512,1024,2048))
    r=x['H_error_ratios_on_bin_doubling']
    assert 3.3<min(r) and max(r)<4.5
    rows=x['mesh_rows']
    assert abs(rows[0]['entropy_error_relative'])>.05
    assert abs(rows[-1]['entropy_error_relative'])<.002
    assert all(abs(z['entropy_identity_error'])<1e-22 for z in rows)
    assert all(abs(z['energy_residual'])<1e-22 for z in rows)
    assert rows[-1]['Qchi1_error_relative']>0 # heat error changes sign on refinement!
    assert rows[-1]['Qchi2_error_relative']<0


def test_invalid_interpolation_remainder_parameters_rejected():
    for h,a,k in ((0.,.3,4.),(.3,-.1,4.),(.3,.3,0.)):
        with pytest.raises(ValueError):
            linear_interpolation_bound(h,a,k)
