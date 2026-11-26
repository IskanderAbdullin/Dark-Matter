"""Independent analytic on-shell support and numerical boundary checks."""
import math
import numpy as np
from kinetics.phi_bose_parent_cutoff_kinematic_ir import (
    daughter_energy_bounds,min_energy_with_parent_cutoff,quadratic_first_cell,
    M1,M2,MPH)


def test_exact_boosted_two_body_energy_extrema_and_global_minimum():
    for p in (0.,.1,1.,2.,8.,16.):
        x=daughter_energy_bounds(p)
        E=math.hypot(p,M1)
        q=x['qstar'];e_star=x['Ephi_star']
        raw_min=(E*e_star-p*q)/M1
        assert abs(raw_min/x['minimum']-1)<1e-8
        for c in np.linspace(-1.,1.,501):
            eph=(E*e_star-p*q*c)/M1
            assert x['minimum']-1e-12<=eph<=x['maximum']+1e-12
    assert daughter_energy_bounds(0.)['minimum']>0.047
    turn=daughter_energy_bounds(2.)['turning_momentum']
    assert 50000.<turn<60000.
    assert abs(daughter_energy_bounds(turn)['minimum']/MPH-1)<1e-8
    assert min_energy_with_parent_cutoff(2.)>min_energy_with_parent_cutoff(8.)>min_energy_with_parent_cutoff(16.)>MPH


def test_n48_n64_zero_first_cell_is_exact_parent_cutoff_support_not_angular_alias():
    t={n:quadratic_first_cell(n,2.) for n in (16,32,40,42,43,48,64)}
    assert 42<t[43]['critical_n']<43
    for n in (16,32,40,42):
        assert not t[n]['first_cell_has_no_possible_on_shell_event']
    for n in (43,48,64):
        assert t[n]['first_cell_has_no_possible_on_shell_event']
    assert abs(t[32]['true_min_daughter_energy_GeV']-0.0011992509375794726)<1e-13
    for pmax in (8.,16.):
        assert not quadratic_first_cell(64,pmax)['first_cell_has_no_possible_on_shell_event']
    # A finite radial/angle sampling can fail despite mathematically allowed
    # support when edge>=minimum. Only edge<minimum proves *no event exists*.
    assert quadratic_first_cell(42,2.)['first_cell_upper_GeV']>min_energy_with_parent_cutoff(2.)
