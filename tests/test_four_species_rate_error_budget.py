"""Targeted cutoff and massive-lepton error gates, NOT physical Γelastic."""
from kinetics.four_species_rate_error_budget import sample


def test_incoming_cutoff_missing_physical_heat_despite_exact_N_E_H():
    a=sample(cut=10.,nmu=6,nc=12,nphi=12)
    b=sample(cut=16.,nmu=6,nc=12,nphi=12)
    assert a['Qchi1_elastic']<0 and b['Qchi1_elastic']<0
    assert a['Qchi2_elastic']>0 and b['Qchi2_elastic']>0
    assert abs(b['Qchi1_elastic']/a['Qchi1_elastic']-1)>.01
    assert a['entropy_elastic']>0 and b['entropy_elastic']>0
    assert max(abs(x['combined_energy_residual']) for x in (a,b))<2e-24
    assert max(x['max_charge_residual'] for x in (a,b))<2e-24


def test_finite_mass_lepton_changes_both_heat_channels_not_yukawa_stoichiometry():
    a=sample(cut=10.,nmu=6,nc=12,nphi=12,mlepton=0.)
    b=sample(cut=10.,nmu=6,nc=12,nphi=12,mlepton=.5)
    assert .85<abs(b['Qchi1_elastic']/a['Qchi1_elastic'])<.95
    assert .85<abs(b['Qchi2_elastic']/a['Qchi2_elastic'])<.95
    assert abs(b['Yukawa_Rphi']/a['Yukawa_Rphi']-1)<2e-12
    assert a['entropy_yukawa']>0 and b['entropy_yukawa']>0
