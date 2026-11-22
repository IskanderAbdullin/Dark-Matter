"""Finite-volume Zprime + distinct nu, antinu event graph: scoped reproducible gates."""
import math
import numpy as np
import pytest
from scipy.special import expit

from kinetics.vector_neutrino_reversible_event_graph import (
    graph,synthetic_spectra,ledger,event_moments,sweep,evolve,
    independent_continuum_weak,
)


def test_exact_thermal_two_chemical_potential_fixed_point():
    g=graph(m=.4,h=.2,emax=4)
    k,E=g["k"],g["E"]
    mu=.1
    fn=expit(mu-k)
    fa=expit(mu-k)
    b=1/np.expm1(E-2*mu)
    after=sweep(g,b,fn,fa,dt=3.)
    for before,computed in zip((b,fn,fa),after):
        assert np.max(np.abs(before-computed))<1e-13
    assert abs(event_moments(g,b,fn,fa)["entropy_production"])<1e-12


def test_eventwise_implicit_nonthermal_strong_evolution_conserves_three_invariants():
    g=graph(m=.4,h=.2,emax=4)
    initial=synthetic_spectra(g)
    start=ledger(g,*initial)
    b,fn,fa,end=evolve(g,*initial,total_t=4.,dt=.25)
    for name,metric in (
        ("reduced_charge","max_abs_reduced_charge_error"),
        ("lepton_asymmetry","max_abs_lepton_asymmetry_error"),
        ("energy","max_abs_energy_error"),
    ):
        assert end[metric]<2e-12*max(1,abs(start[name]))
    assert end["entropy_gain"]>0 and end["min_step_entropy_gain"]>0
    assert np.min(b)>=0
    assert np.min(fn)>=0 and np.min(fa)>=0
    assert np.max(fn)<=1 and np.max(fa)<=1
    assert not np.array_equal(fn,fa)


def test_from_empty_vector_has_quantum_inverse_and_preserves_bounds():
    g=graph(.4,.2,4)
    b=np.zeros_like(g["E"])
    fn=expit(-g["k"])
    fa=expit(-g["k"])
    prior=ledger(g,b,fn,fa)
    out=sweep(g,b,fn,fa,dt=3.)
    after=ledger(g,*out)
    assert np.min(out[0])>0
    assert 0<=np.min(out[1])<=np.max(out[1])<=1
    assert 0<=np.min(out[2])<=np.max(out[2])<=1
    assert after["entropy"]>prior["entropy"]
    for metric in ("energy","reduced_charge","lepton_asymmetry"):
        assert math.isclose(after[metric],prior[metric],rel_tol=1e-12,abs_tol=1e-12)


def test_resolution_alias_cannot_claim_continuum_at_coarse_h():
    m=.4;emax=4
    ref=independent_continuum_weak(m,emax,nr=144,nangle=144)
    coarse=event_moments(graph(m,.4,emax),
                         *synthetic_spectra(graph(m,.4,emax)))
    fine_graph=graph(m,.05,emax)
    fine=event_moments(fine_graph,*synthetic_spectra(fine_graph))
    assert abs(coarse["dNV"]/ref["dNV"]-1)>.1
    assert abs(coarse["dEV"]/ref["dEV"]-1)>.1
    # Finite-volume charge/energy exactness does NOT imply full continuum
    # quadrature accuracy; the fine benchmark still has ~0.4% errors.
    for metric in ("dNV","dEV","entropy_production"):
        assert abs(fine[metric]/ref[metric]-1)<.01
    assert abs(fine["dEV"]/ref["dEV"]-1)>1e-3


def test_time_discretization_has_real_bias_despite_exact_event_conservation():
    g=graph(.4,.2,4)
    init=synthetic_spectra(g)
    trials={}
    for dt in (1.,.5,.25,.125):
        trials[dt]=evolve(g,*init,total_t=4.,dt=dt)[0]
    ref=trials[.125]
    err1=np.linalg.norm(trials[1.]-ref)/np.linalg.norm(ref)
    err05=np.linalg.norm(trials[.5]-ref)/np.linalg.norm(ref)
    err025=np.linalg.norm(trials[.25]-ref)/np.linalg.norm(ref)
    assert err1>err05>err025>1e-3


def test_reject_bad_energy_lattice_and_pauli():
    with pytest.raises(ValueError):
        graph(m=.3,h=.2,emax=4)
    g=graph()
    b,fn,fa=synthetic_spectra(g)
    fn[0]=1.01
    with pytest.raises(ValueError):
        sweep(g,b,fn,fa,dt=.1)
