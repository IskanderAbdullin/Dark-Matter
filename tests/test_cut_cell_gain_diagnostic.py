"""Unbiased source-only cut-cell quadrature around hyperbolic on-shell edge."""
from __future__ import annotations
import math
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import dblquad
from scipy.special import expit

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from cut_cell_gain_diagnostic import integrate_cut_cells,evaluate_cells,run
from uv_source_decomposition import exact_rectangle_source


def test_cut_cell_integral_reproduces_independent_direct_double_quadrature():
    edges=np.array([0.,.03,.06,.10,.20])
    m=.1
    source=integrate_cut_cells(edges,mass=m,quadrature_n=32)
    c=m*m/4
    normal=4*m/(2*math.pi**2)
    for i,j in ((1,2),(1,3),(2,2),(2,3),(3,3)):
        xl,xh=edges[i:i+2]
        yl,yh=edges[j:j+2]
        lower=max(xl,c/yh)
        if lower>=xh:
            assert source[i,j]==0
        else:
            direct=dblquad(
                lambda y,x:float(expit(-(x+y))),
                lower,xh,lambda x:max(yl,c/x),lambda x:yh,
                epsabs=1e-12,epsrel=1e-11
            )[0]
            assert abs(source[i,j]-direct)<max(1e-12,1e-8*direct)
    assert abs(normal*source.sum()/exact_rectangle_source(edges[-1],m)-1)<1e-9


def test_full_mesh_source_components_close_exactly():
    a=evaluate_cells(max_index=5000)
    assert a["topology_active_collision_invariants"]==3
    assert a["daughter_modes"]==199
    assert abs(a["relative_cut_cell_quadrature_error"])<2e-6
    assert abs(a["error_decomposition_residual"])<1e-8
    assert a["kinematically_forbidden_centroid_cell_integral"]>=0
    assert a["sampled_kinetic_vacuum_source"]>0


def test_more_gauss_nodes_converge_on_same_cutoff():
    a=evaluate_cells(max_index=5000,quadrature_n=16)
    b=evaluate_cells(max_index=5000,quadrature_n=32)
    assert abs(a["all_cut_cell_integral"]/b["all_cut_cell_integral"]-1)<1e-6
    assert abs(b["relative_cut_cell_quadrature_error"])<2e-6


def test_diagnostic_does_not_modify_collisions_or_invent_class_status():
    a=run(quick=True)
    assert a["no_collision_operator_change"]
    assert len(a["cases"])==1
    assert a["relative_change_gauss16_to_gauss32"]<1e-6
    assert "full nonlinear" in a["remaining"]
    assert "Continuum thermal VACUUM" in a["cases"][0]["warning"]


@pytest.mark.parametrize("mass,nodes",[(0.,16),(.1,1)])
def test_invalid_quadrature_input_rejected(mass,nodes):
    with pytest.raises(ValueError):
        integrate_cut_cells([0,.1,.2],mass=mass,quadrature_n=nodes)
