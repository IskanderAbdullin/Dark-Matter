"""Nonthermal 5D Pauli χℓ weak energy/entropy vs existing independent FD 5D."""
import math
import numpy as np
from scipy.special import expit
from kinetics.nonthermal_chi_lepton_5d_weak import event,weak_moments
from kinetics.finite_mediator_5d_adaptive_reference import full_heat


def test_numerical_full_5d_fd_subset_matches_preexisting_thermal_reference():
    for m in (0.,.5):
        f=lambda E:expit(-np.asarray(E)/.9)
        g=lambda E:expit(-np.asarray(E))
        new=weak_moments(f,g,M=.1,m=m,Te=1.,G=.001,mV=.1,
                         nr=4,nmu=4,nphi=12,nang=12,cut_kinetic_T=12.)
        original=full_heat(M=.1,m=m,Te=1.,Tchi=.9,G=.001,mV=.1,
                           nr=4,nmu=4,nphi=12,polar="cosine",nang=12,
                           cut_kinetic_T=12.)
        assert math.isclose(new["Qchi"],original,rel_tol=2e-11)
        assert new["Qchi"]>0 and new["Qlepton"]==-new["Qchi"]
        assert new["entropy_production"]>0


def test_common_temperature_with_independent_chemical_potentials_is_exact_fixed_point():
    for m in (0.,.5):
        f=lambda E:expit(.2-np.asarray(E))
        g=lambda E:expit(-.1-np.asarray(E))
        x=weak_moments(f,g,m=m)
        noneq=weak_moments(lambda E:expit(-np.asarray(E)/.9),
                           lambda E:expit(-np.asarray(E)),m=m)
        assert abs(x["Qchi"])<1e-13*abs(noneq["Qchi"])
        assert abs(x["entropy_production"])<1e-12*noneq["entropy_production"]


def test_nonthermal_finite_vector_Pauli_energy_and_entropy_on_every_event():
    f=lambda E:expit(-np.asarray(E)/.9+.15*np.exp(-np.asarray(E)/.9))
    g=lambda E:expit(-np.asarray(E)+.05*np.exp(-np.asarray(E)))
    angle=2*np.pi*np.arange(16)/16
    for M,m in ((.1,0.),(.02,.5)):
        for c in (-.85,.1,.9):
            ev=event(M,m,.9,1.2,-.35,.001,.1,c,angle,f,g)
            assert np.max(abs(ev["dEchi"]+ev["dElepton"]))<3e-14
            assert np.all(ev["entropy"]>=-1e-25)
            assert ev["amplitude"]>=0
    x=weak_moments(f,g,m=0.,nphi=12)
    assert x["Qchi"]>0 and x["entropy_production"]>0
    assert math.isclose(x["Qlepton"],-x["Qchi"],abs_tol=0)


def test_interior_fd_validation_fail_closed():
    bad=lambda E:np.ones_like(np.asarray(E))
    normal=lambda E:expit(-np.asarray(E))
    import pytest
    with pytest.raises(ValueError):
        event(.1,0,.9,1.2,-.35,.001,.1,.1,np.array([0.,1.]),bad,normal)
