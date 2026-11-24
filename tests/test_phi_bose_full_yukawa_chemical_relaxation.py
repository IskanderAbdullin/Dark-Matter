"""Physically normalized near-equilibrium FIVE Dirac/anti/phi Yukawa chemistry.

Restricted fixed-T FD/BE fugacity manifold is NOT whole phase-space
thermalization, FRW or 86eV relic Ω; the true conserved vector U1 is exact.
"""
import math
import numpy as np
from kinetics.phi_bose_full_yukawa_chemical_relaxation import (
    susceptibility,equilibrium_event_rate_weight,
    equilibrium_chemical_relaxation,Q_U1,
)
from kinetics.phi_bose_five_dirac_full_yukawa_implicit import (
    FullDiracPhysicalMoments
)
from kinetics.phi_bose_yukawa_charge_audit import stoichiometry


def _model(P=2):
    return FullDiracPhysicalMoments(P=P,Ppair=P,N1=3,N2=3,Nphi=8,
                      nrad=4,nangle=8,pair_nrad=4,
                      pair_nrelative=8,pair_ncm=32,pair_naz=4)


def test_independent_thermal_susceptibility_massless_FD_Bose_limits():
    # One fermion species including two spin states vs one real Bose state.
    for T in (.1,1.,2.):
        fd=susceptibility(T,0.,2,"fermion")
        be=susceptibility(T,0.,1,"boson")
        assert math.isclose(fd,T*T/6,rel_tol=5e-11)
        assert math.isclose(be,T*T/6,rel_tol=5e-11)
        assert susceptibility(T,.02,2,"fermion")<fd
        assert susceptibility(T,86e-9,1,"boson")<be


def test_chemical_linear_response_is_nonnegative_and_one_true_U1_zero_mode():
    m=_model()
    q=equilibrium_chemical_relaxation(m,T=1.,gstar=61.75)
    evals=np.asarray(q["chemical_relaxation_eigenvalues_GeV"])
    assert evals.shape==(5,)
    assert evals[-1]>0
    assert abs(evals[0])<1e-10*evals[-1]
    assert all(e>1e-9*evals[-1] for e in evals[1:])
    assert q["true_vector_U1_symmetric_null_residual"]<1e-11
    s=stoichiometry(True)
    assert np.array_equal(s@Q_U1,np.zeros(4))
    chi=np.asarray(q["susceptibility_GeV2"])
    gamma=np.asarray(q["channel_equilibrium_F_integrated_GeV4"])
    L=(s.T@np.diag(gamma)@s)/np.sqrt(np.outer(chi,chi))
    for v in ([1,2,-2,4,-1],[-1,0,1,0,0],[0,0,0,0,1]):
        assert np.asarray(v)@L@np.asarray(v)>=-1e-24


def test_equilibrium_microevent_weights_and_finite_difference_chemical_susceptibility():
    m=_model()
    T=1.
    a=equilibrium_chemical_relaxation(m,T)
    gamma=np.asarray(a["channel_equilibrium_F_integrated_GeV4"])
    s=stoichiometry(True)
    direction=np.array([.02,-.01,.01,-.025,-.005])
    predicted=-(gamma/T)*(s@direction)
    measured=[]
    for epsilon in (1e-3,5e-4):
        mu=direction*epsilon
        e=m.initial_etas(T=(T,)*5,mu=mu)
        source=m.source(e)
        observed=np.array([c["event_rate"] for c in source["channels"]])
        assert np.all(np.isfinite(observed))
        assert np.linalg.norm(observed)>0
        assert np.max(abs(observed/epsilon-predicted)/
                      np.maximum(abs(predicted),1e-24))<5e-3
        measured.append(observed/epsilon)
    assert np.linalg.norm(measured[1]-predicted)<=(
           np.linalg.norm(measured[0]-predicted)*.7)


def test_scalar_pair_channels_repair_old_accidental_number_charges_and_H_eigenmode():
    m=_model()
    eq=equilibrium_chemical_relaxation(m,T=1.)
    gamma=np.asarray(eq["channel_equilibrium_F_integrated_GeV4"])
    assert all(g>0 for g in gamma)
    assert gamma[0]==gamma[1]  # μU=0 exact conjugation
    assert gamma[2]>0 and gamma[3]>0
    # The positive pair weights remove two accidental decay-only
    # number charges in the full FIVE-species stoichiometric matrix.
    s=stoichiometry(True)
    assert np.linalg.matrix_rank(s)==4
    assert s.shape[1]-np.linalg.matrix_rank(s)==1
    assert eq["chemical_relaxation_over_H"][1]>0


def test_cutoff_dependence_is_measured_not_assumed_converged():
    a=equilibrium_chemical_relaxation(_model(2.),T=1.)
    b=equilibrium_chemical_relaxation(_model(4.),T=1.)
    ga=np.asarray(a["channel_equilibrium_F_integrated_GeV4"])
    gb=np.asarray(b["channel_equilibrium_F_integrated_GeV4"])
    assert np.all(ga>0) and np.all(gb>0)
    # Restricting high parent/pair momenta changes absolute rates;
    # graph conservation/H cannot validate a truncated collision integral.
    assert np.max(abs(gb/ga-1))>1e-3
    assert a["true_vector_U1_symmetric_null_residual"]<1e-11
    assert b["true_vector_U1_symmetric_null_residual"]<1e-11
