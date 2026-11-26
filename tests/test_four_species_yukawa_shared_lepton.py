"""Scoped common-bath χ1/χ2/real-φ/lepton Pauli/Bose kinetics, not 86eV DM."""
import math
import numpy as np
import pytest
from scipy.special import expit
from kinetics.four_species_yukawa_shared_lepton import (
    configuration,initial,rhs,ledger,evolve
)


def test_four_species_combined_source_has_independent_quantum_entropy_and_energy():
    c=configuration(bins=16)
    x=initial(c,T1=1.1,T2=1.,Tlep=1.,Tphi=1.)
    ds,obs=rhs(c,x)
    assert obs['elastic_entropy']>0 and obs['yukawa_entropy']>0
    assert obs['entropy_event']>0 and obs['entropy_node']>0
    assert obs['entropy_identity_abs_error']<3e-24
    assert abs(obs['energy_residual'])<3e-24
    assert max(abs(q) for q in obs['charge_residuals'])<3e-24
    assert obs['Yukawa_Rphi']>0 and obs['elastic_Qlep']>0
    assert sum(obs['elastic_Qchi'])<0
    assert abs(obs['Yukawa_Qtotal'])<3e-24
    e=c['elastic'];y=c['yukawa']
    assert np.array_equal(y.w1,e['chi_w'][0])
    assert np.array_equal(y.w2,e['chi_w'][1])


@pytest.mark.parametrize('mlepton',(0.,.5))
def test_all_four_quantum_species_common_T_with_independent_chemical_mu_stationary(mlepton):
    c=configuration(bins=12,mlepton=mlepton)
    e=c['elastic'];y=c['yukawa'];T=.8
    state=(expit((.03-y.e1)/T),expit((.05-y.e2)/T),
           expit((.07-e['lg'])/T),1/np.expm1((y.ephi+.02)/T))
    ds,z=rhs(c,state)
    assert max(np.max(abs(q)) for q in ds)<1e-19
    assert abs(z['entropy_event'])<1e-32


def test_nonlinear_four_species_time_preserves_three_charges_and_entropy():
    c=configuration(bins=16)
    x=initial(c,T1=1.1,T2=1.,Tlep=1.)
    z=evolve(c,x,end_time=2e7,initial_dt=1e7,atol=2e-5)
    assert z['accepted']>=4 and z['rejected']>=1
    assert z['Nphi_created']>0 and z['N2_created']>0
    assert math.isclose(z['Nphi_created'],z['N2_created'],abs_tol=1e-13)
    assert max(abs(z[k]) for k in ('dNdark','dN1phi','dNlep','dEnergy'))<2e-12
    assert z['dElepton']>0
    assert z['min_step_entropy_gain']>0 and z['dEntropy']>0
    assert 0<z['Fermi_min']<=z['Fermi_max']<1
    assert z['Bose_min']>0 and z['max_local_spectral_error']<2e-5


def test_disabled_elastic_coupling_yields_no_lepton_heat_but_yukawa_survives():
    c=configuration(bins=16,elastic_G=0.)
    x=initial(c,T1=1.1)
    ds,z=rhs(c,x)
    assert abs(z['elastic_entropy'])<1e-26
    assert abs(z['elastic_Qlep'])<1e-26
    assert z['Yukawa_Rphi']>0
    out=evolve(c,x,end_time=1e7,initial_dt=1e7,atol=2e-5)
    assert abs(out['dElepton'])<1e-15
    assert out['Nphi_created']>0
    assert out['dEntropy']>0


def test_yukawa_daughters_must_fit_shared_energy_grid_and_bose_zero_fail_closed():
    with pytest.raises(ValueError):
        configuration(bins=16,elastic_cut=1.,yukawa_pmax=12.)
    c=configuration(bins=16)
    x=list(initial(c))
    x[3]=np.zeros_like(x[3])
    with pytest.raises(ValueError):
        rhs(c,tuple(x))
