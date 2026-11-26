"""Microreversible F/R chemical fixed point can coexist with wrong thermal Nφ."""
from kinetics.phi_bose_false_equilibrium_gate import compare


def test_chemical_fixed_point_does_not_certify_scalar_number_86ev():
    for n in (16,64,256):
        x=compare(n,86e-9)
        assert x['absolute_yukawa_rate']<1e-20
        assert x['absolute_entropy_production']<1e-24
        assert x['absolute_energy_balance']<1e-20
        assert x['max_absolute_charge_defect']<1e-20
        assert x['false_full_thermal_number_ratio']>10
        assert x['false_first_cell_ratio']>1e5


def test_mass_and_resolution_change_false_first_cell_without_BEC():
    x=compare(64,.001)
    y=compare(64,86e-9)
    assert y['false_full_thermal_number_ratio']>1000
    assert 1.1 < x['false_full_thermal_number_ratio']<1.2
    assert y['false_first_cell_ratio'] > 1e4*x['false_first_cell_ratio']
