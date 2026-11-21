"""Regridding feasibility: independent Pauli-bounded energy interval and repairs."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pytest
from scipy.optimize import linprog
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from pauli_energy_bounds import (
    shell_capacities_and_energy, feasible_energy_interval,
    repair_energy_by_bounded_transfer, run
)
from parent_q_remap import number, energy, rebin


def test_pauli_extrema_agree_with_independent_linear_programming():
    q=np.array([0.,.03,.1,.25,.8,1.5,3.])
    cap,E=shell_capacities_and_energy(q,.2)
    particles=cap.sum()*.37
    interval=feasible_energy_interval(q,.2,particles)
    # Independent reference: linear programming in per-shell particle number.
    low=linprog(E,A_eq=np.ones((1,len(cap))),b_eq=[particles],
                bounds=[(0,float(c)) for c in cap],method="highs")
    high=linprog(-E,A_eq=np.ones((1,len(cap))),b_eq=[particles],
                 bounds=[(0,float(c)) for c in cap],method="highs")
    assert low.success and high.success
    assert np.isclose(interval["minimum_energy"],low.fun,rtol=1e-10,atol=1e-12)
    assert np.isclose(interval["maximum_energy"],-high.fun,rtol=1e-10,atol=1e-12)


def test_corrected_two_moment_remap_is_bounded_but_changes_shape():
    old=np.r_[0.,np.geomspace(.01,5.,35)]
    new=np.r_[0.,np.geomspace(.012,5.,43)]
    f=1/(1+np.exp(np.hypot(.5*(old[1:]+old[:-1]),.1)))
    x=.3
    result=repair_energy_by_bounded_transfer(old,new,f,x)
    final=result["corrected_occupations"]
    assert np.all((final>=0)&(final<=1))
    assert abs(number(new,final)-number(old,f))<1e-11
    assert abs(energy(new,final,x,gauss_n=24)-energy(old,f,x,gauss_n=24))<1e-9
    assert abs(result["number_error"])<1e-11
    assert abs(result["relative_energy_error"])<1e-9
    assert np.isfinite(result["numeric_entropy_change"])
    assert "NOT PHYSICAL" in result["classification"]


def test_energy_infeasible_if_new_grid_one_cell_fixes_mean_energy():
    old=[0.,.01,.02,.2]
    new=[0.,.2]
    f=[1.,0.,0.]
    N=number(old,f)
    interval=feasible_energy_interval(new,.1,N)
    E=energy(old,f,.1,gauss_n=24)
    assert E<interval["minimum_energy"]
    with pytest.raises(ValueError,match="infeasible"):
        repair_energy_by_bounded_transfer(old,new,f,.1)


def test_empty_population_energy_interval_is_zero():
    a=feasible_energy_interval([0.,1.,2.],.1,0.)
    assert a["minimum_energy"]==0
    assert a["maximum_energy"]==0


def test_filling_all_states_energies_collapse_to_one_value():
    q=[0.,.1,.5,2.]
    cap,_=shell_capacities_and_energy(q,.1)
    interval=feasible_energy_interval(q,.1,float(cap.sum()))
    assert abs(interval["minimum_energy"]-interval["maximum_energy"])<1e-10


def test_invalid_pauli_demand_is_rejected():
    cap,_=shell_capacities_and_energy([0,1,2],.5)
    with pytest.raises(ValueError):
        feasible_energy_interval([0,1,2],.5,float(cap.sum()*1.01))
    with pytest.raises(ValueError):
        feasible_energy_interval([0,1,2],.5,-1)


def test_cli_labels_energy_correction_as_diagnostic_not_cosmology():
    result=run()
    assert "NOT PHYSICAL" in result["classification"]
    assert abs(result["number_error"])<1e-10
    assert abs(result["relative_energy_error"])<1e-9
