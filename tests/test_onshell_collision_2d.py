"""Independent physical on-shell 2D collision-integral regression tests."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
from scipy.integrate import quad
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from onshell_collision_2d import OnShellCollision2D, allowed_rectangle_areas, vacuum_reference
from one_to_two_fixed_bath import FlatThermalDecay


def test_analytic_continuum_vacuum_source_matches_original_decay_quadrature():
    original=FlatThermalDecay(nq=50,np1=140,nmu=72)
    rate=float(np.dot(original.rateprefactor,original.parent_fermi))
    exact=vacuum_reference()
    assert abs(rate/exact-1)<1e-8


def test_hyperbola_area_is_integrated_analytically():
    edges=np.geomspace(1e-5,32,60)
    areas=allowed_rectangle_areas(edges)
    c=.1**2/4
    rectangle=np.diff(edges)[:,None]*np.diff(edges)[None,:]
    mask=(areas>0)&(areas<rectangle*(1-1e-5))
    assert mask.sum()>10
    ij=np.argwhere(mask)
    for i,j in ij[::max(1,len(ij)//8)]:
        xl,xh=edges[i:i+2]
        yl,yh=edges[j:j+2]
        f=lambda x:max(0.,yh-max(yl,c/x))
        ref=quad(f,xl,xh,epsabs=1e-13,points=[x for x in (c/yh,c/yl) if xl<x<xh])[0]
        assert abs(ref-areas[i,j])<max(1e-12,1e-6*ref)


@pytest.mark.parametrize("n",[50,100,150])
def test_all_active_collision_channels_are_on_shell(n):
    m=OnShellCollision2D(n)
    i,j=m.forward.nonzero()
    assert i.size>0
    assert np.min(4*m.e[i]*m.e[j]-.1**2)>=-1e-13
    assert m.max_allowed_pair_mass_deficit==0.
    assert 0<m.vacuum_rate<m.vacuum_exact


def test_vacuum_cell_integration_converges_toward_exact_1d_reference():
    coarse=OnShellCollision2D(50)
    medium=OnShellCollision2D(100)
    fine=OnShellCollision2D(220)
    assert abs(fine.vacuum_relative_error)<abs(medium.vacuum_relative_error)
    assert abs(medium.vacuum_relative_error)<abs(coarse.vacuum_relative_error)
    assert abs(fine.vacuum_relative_error)<.01


@pytest.mark.parametrize("mu",[0.,-.3,-2.])
def test_on_shell_collision_discrete_thermal_stationarity(mu):
    m=OnShellCollision2D(80)
    y=m.equilibrium(mu)
    d=m.rhs(0,y)
    assert np.max(abs(d[:m.nq]*m.vol_phi))<1e-11
    assert np.max(abs(d[m.nq:]*m.vol_chi))<1e-11


def test_exact_daughter_invariant_bounds_and_discrete_h_theorem():
    m=OnShellCollision2D(100)
    checks=m.checks()
    assert abs(checks["daughter_number_rate_mismatch"])<1e-11
    assert checks["free_energy_rate"]<0
    assert checks["entropy_production"]>0
    assert abs(checks["dF_plus_sigma"])<1e-10


def test_on_shell_burst_number_and_mean_momentum_have_partial_grid_convergence():
    a=OnShellCollision2D(100).evolve()
    b=OnShellCollision2D(150).evolve()
    assert abs(a["daughter_difference_over_T3"])<2e-7
    assert abs(b["daughter_difference_over_T3"])<2e-7
    assert b["nphi_over_T3"]>a["nphi_over_T3"]>0.3
    assert abs(b["mean_phi_q"]/a["mean_phi_q"]-1)<.03
    assert b["max_fermion_occupation"]<=1.000001


def test_on_shell_pair_flux_energy_and_number_ledgers():
    m=OnShellCollision2D(100)
    # Nonthermal but physically admissible Bose/fermion occupations.
    y=np.r_[0.1*np.exp(-m.e/.2),0.25*np.exp(-m.e)]
    b=m.moment_budget(y)
    assert abs(b["daughter_number_difference_rate"])<1e-11
    assert abs(b["phi_number_rate"]-b["reaction_number_rate"])<1e-11
    assert abs(b["chi2_number_rate"]-b["reaction_number_rate"])<1e-11
    assert abs(b["on_shell_energy_ledger_residual"])<1e-11
    assert abs(b["daughter_total_energy_rate"]+
               b["on_shell_parent_reservoir_energy_loss"])<1e-11
