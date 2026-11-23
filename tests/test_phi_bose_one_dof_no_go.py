"""Focused one-DOF physical Bose number/energy/fixed-point falsification."""
import math
import numpy as np
from kinetics.phi_bose_one_dof_no_go import cells,audit


def test_one_eta_monotone_number_and_locked_energy():
    e,w,n,u,eta,ebar=cells(32,math.hypot(2.,.1)+.2)
    occupancy=1/np.expm1(eta)
    assert np.allclose(n,w*occupancy,rtol=3e-15,atol=1e-16)
    assert np.all(-w*occupancy*(1+occupancy)<0)
    assert np.all((e[:-1] < ebar) & (ebar < e[1:]))
    assert abs(float(np.sum(n))-.056268555792202)<1e-11
    assert np.sum(ebar*n-u)>0


def test_false_common_temperature_equilibrium_and_H():
    r=audit(32)
    assert abs(r['exact_thermal_scalar_rate'])<1e-22
    assert r['spurious_equilibrium_scalar_rate']>1e-14
    assert r['spurious_entropy_rate']>1e-15
    assert r['scalar_thermal_U_relative_error']>5e-4
    assert r['energy_closure_residual']>1e-15
    assert r['residual_over_absolute_energy_throughput']>.02


def test_three_mode_refinements_reduce_error_not_exactly_remove():
    x=[audit(n) for n in (16,32,64)]
    assert x[0]['scalar_thermal_U_relative_error']>x[1]['scalar_thermal_U_relative_error']>x[2]['scalar_thermal_U_relative_error']>0
    assert x[0]['spurious_entropy_rate']>x[1]['spurious_entropy_rate']>x[2]['spurious_entropy_rate']>0
    assert all(v['class_ready'] is False for v in x)
