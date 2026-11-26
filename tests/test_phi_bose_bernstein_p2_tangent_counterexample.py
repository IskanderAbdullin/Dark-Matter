"""Physical quadratic Bose tangent can turn negative inside an on-shell cell.

Independently checks scalar-energy-first nonthermal-mode normalization
against parent-first three-field common event source. A negative projected
inverse-eta derivative does NOT mean the physical occupation becomes negative
or that a condensate forms.
"""
import math

import numpy as np

from kinetics.phi_bose_bernstein_p2_tangent_counterexample import (
    first_active_cell, independent_bernstein_p2_tangent,
    event_order_bernstein_tangent, true_eta_dot, projected_minimum,
)


def test_quadratic_bose_projection_false_negative_inside_real_supported_cell():
    for P,N in ((8.,8),(16.,8),(16.,12)):
        x=independent_bernstein_p2_tangent(P,N,64,64)
        mini=x["projected_minimum"]
        assert x["min_eigenvalue_G"]>0
        assert mini["value"]<-3.5e-10,(P,N,mini)
        assert 0.6<mini["x"]<0.7
        assert true_eta_dot(mini["E_GeV"],P,64)>0
        L,U=x["energy_interval"]
        for k in np.linspace(L+(U-L)*1e-4,U-(U-L)*1e-4,11):
            assert true_eta_dot(float(k),P,64)>0
        assert x["d_eta_coefficients"][0]>0
        assert x["d_eta_coefficients"][1]<0
        assert x["d_eta_coefficients"][2]>0
        assert x["physical_firstcell_Rphi"]<0
        assert x["physical_firstcell_Qphi"]<0


def test_P8_N8_false_negative_and_refinement_reproducible():
    values=[]
    for order in (32,48,64,80):
        result=independent_bernstein_p2_tangent(8.,8,order,order)
        values.append(result["projected_minimum"]["value"])
    assert all(x<0 for x in values)
    assert math.isclose(values[-1],-5.159774281937157e-10,rel_tol=6e-4)
    assert abs(values[-1]-values[-2])/abs(values[-1])<3e-3
    assert math.isclose(
        true_eta_dot(.08397890951664282,8.,80),
        9.68624932417159e-10,rel_tol=1e-5)


def test_negative_middle_Bernstein_coefficient_is_not_by_itself_a_negative_polynomial():
    for P,N in ((2.,8),(8.,12),(8.,16),(16.,16)):
        result=independent_bernstein_p2_tangent(P,N,64,64)
        coeff=np.array(result["d_eta_coefficients"])
        assert coeff[1]<0
        assert result["projected_minimum"]["value"]>0
        L,U=result["energy_interval"]
        assert math.isclose(
            result["projected_minimum"]["value"],
            min(coeff[0],coeff[2],
                (lambda x: coeff[0]*(1-x)**2+
                  2*coeff[1]*x*(1-x)+coeff[2]*x*x)(
                      result["projected_minimum"]["x"])),
            rel_tol=2e-12)


def test_higher_resolution_removes_specified_P8_N8_inside_negative():
    for P,old_N,new_N in ((8.,8,16),(16.,8,16),(16.,12,24)):
        coarse=independent_bernstein_p2_tangent(P,old_N,64,64)
        fine=independent_bernstein_p2_tangent(P,new_N,64,64)
        assert coarse["projected_minimum"]["value"]<0
        assert fine["projected_minimum"]["value"]>0


def test_independent_parent_first_same_physical_energy_tangent_sign():
    for P,N in ((2.,8),(8.,8),(8.,16),(16.,16)):
        independent=independent_bernstein_p2_tangent(P,N,64,64)
        common=event_order_bernstein_tangent(P,N,8,16,32)
        a=np.array(independent["d_eta_coefficients"])
        b=np.array(common["d_eta_coefficients"])
        assert np.max(abs((a-b)/np.maximum(abs(a),3e-10)))<2e-2,(P,N,a,b)
        assert math.copysign(1.,independent["projected_minimum"]["value"])==math.copysign(
            1.,common["projected_minimum"]["value"])
