"""Quantum exact shell graph: physical kinematics, nullspace, H, bounds, time."""
import numpy as np
from kinetics.phi_bose_exact_shell_graph import graph, invariant_audit, equilibrium, implicit_event, network_run


def test_kinematically_allowed_exact_energy_and_graph_nullity():
    for n in (2,4,8,16):
        energies,events=graph(n)
        for i,j,k,c in events:
            assert abs(energies[0][i]-energies[1][j]-energies[2][k])<1e-12
            assert abs(c)<=1
        x=invariant_audit(n)
        assert x['left_nullity']==3
        assert x['stoichiometric_rank']==3*(n+1)-3
        assert x['invariant_residual']<1e-12


def test_true_thermal_common_temperature_chemical_fixed_point():
    for n in (2,8,16):
        x=equilibrium(n)
        assert x['max_relative_difference']<1e-13


def test_positive_bose_pauli_implicit_extent_H_and_exact_shell_energy():
    energies,_=graph(4)
    x=implicit_event(.5,.2,.7,dt=3,rate=.7,capacities=(1.1,.8,1.2))
    a,d,b=x['new']
    assert 0<a<1 and 0<d<1 and b>0
    assert abs(x['implicit_residual'])<2e-15
    assert x['entropy_gain']>=x['entropy_lower_bound']-2e-15
    assert x['entropy_gain']>0
    u=implicit_event(.5,.2,.7,dt=3,rate=.7)
    assert u['entropy_gain']>=u['entropy_lower_bound']-2e-15
    assert abs(u['extent']*(-energies[0][2]+energies[1][2]+energies[2][0]))<1e-13


def test_full_resonance_graph_time_preserves_ledger_and_bounds():
    runs=[network_run(8,nsteps=n) for n in (1,2,4,8)]
    for r in runs:
        assert max(abs(x) for x in r['invariants_residual'])<1e-13
        assert r['entropy_gain']>0
        assert 0<r['bounds']['f1_min']<=r['bounds']['f1_max']<1
        assert 0<r['bounds']['f2_min']<=r['bounds']['f2_max']<1
        assert r['bounds']['phi_min']>0
    assert abs(runs[0]['delta_Nphi']-runs[-1]['delta_Nphi'])>1e-6
