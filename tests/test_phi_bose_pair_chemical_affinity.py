"""Nonzero chemical potential fixed point of restricted decay is NOT full real-phi Yukawa equilibrium."""
import math
import numpy as np
import pytest
from kinetics.phi_bose_pair_chemical_affinity import (
    chemical_affinities,charge_symmetric_potentials,
    paired_fermion_bose_weights,event_entropy,
)
from kinetics.phi_bose_yukawa_charge_audit import stoichiometry


def test_full_network_chemical_equilibrium_has_only_real_vector_U1():
    S=stoichiometry(True)
    mu=charge_symmetric_potentials(.3,.3,0.)
    assert np.array_equal(S@mu,np.zeros(4))
    assert np.array_equal(chemical_affinities(mu),[0.,0.,0.,0.])
    assert np.linalg.matrix_rank(S)==4
    # Exact single continuous null charge for full 4-reaction/5-species network.
    assert S.shape[1]-np.linalg.matrix_rank(S)==1
    # With antiparticle mu=-particle mu, the pair channels require
    # 2 mu_phi=0, and the two decay channels then require mu1=mu2.
    mu2=charge_symmetric_potentials(-.15,.05,-.2)
    a=np.asarray(chemical_affinities(mu2))
    assert np.max(abs(a-[0.,.4,.4,.4]))<1e-14


def test_real_scalar_pair_microevent_detailed_balance_exact_chemical_affinity():
    T=1.
    for mu_phi in (-.2,-.02,-.001,0.):
        for e2,ea,k1,k2 in ((.1,.1,.09,.11),
                             (.6,.3,.4,.5),
                             (.3,.7,.3,.7)):
            F,R=paired_fermion_bose_weights(
                e2,ea,k1,k2,.05,-.05,mu_phi,T)
            assert F>0 and R>0
            assert math.isclose(F/R,math.exp(-2*mu_phi/T),
                                rel_tol=2e-13)
            assert event_entropy(F,R)>=-1e-14
            if mu_phi==0.:
                assert math.isclose(F,R,rel_tol=2e-13)


def test_nonzero_scalar_mu_restricted_particle_fixed_point_has_nonzero_pair_flux():
    F,R=paired_fermion_bose_weights(.1,.1,.09,.11,.05,-.05,-.2)
    assert F>R and event_entropy(F,R)>0
    assert math.isclose(F/R,math.exp(.4),rel_tol=2e-13)
    # Previous particle-only chi1->chi2+phi chemical equilibrium (A=0)
    # also has a NONstationary conjugate antiparticle decay.
    assert abs(chemical_affinities(
        charge_symmetric_potentials(-.15,.05,-.2))[1]-.4)<1e-13


def test_chemical_occupancy_validation_and_microphysical_scope():
    with pytest.raises(ValueError):
        paired_fermion_bose_weights(.1,.2,.1,.1,.05,-.05,-.2)
    with pytest.raises(ValueError):
        paired_fermion_bose_weights(.1,.2,.1,.2,.05,-.05,.2)
    with pytest.raises(ValueError):
        chemical_affinities([0.,0.])
