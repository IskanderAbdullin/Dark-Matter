"""Quantum bounded implicit two-level reduced finite-bath kinetic tests."""
from __future__ import annotations
import random,sys
from pathlib import Path
import pytest
from scipy.special import expit
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from fermion_finite_bath_pairs import evolve_pair,evolve_network


def test_exact_FD_and_finite_bath_stationary_for_arbitrary_large_dt():
    A=10.;T=1.5;Ei=.2;Ej=1.3;mu=.25
    fi=expit((mu-Ei)/T);fj=expit((mu-Ej)/T)
    for dt in (1e-4,1,1e5):
        r=evolve_pair(Ei,Ej,fi,fj,3.,4.,A*T**4,A,2.5,dt)
        assert abs(r["dN"])<1e-12
        assert abs(r["entropy_delta"])<1e-12
        assert abs(r["energy_residual"])<1e-12


def test_large_timesteps_never_break_Pauli_charge_heat_or_entropy():
    for dt in (.01,1.,100.,1e5):
        r=evolve_pair(.1,1.,.9,.1,3.,5.,50.,10.,2.,dt)
        assert 0<=r["fi"]<=1 and 0<=r["fj"]<=1 and r["ESM"]>=0
        assert abs(r["particle_residual"])<1e-12
        assert abs(r["energy_residual"])<1e-12
        assert r["entropy_delta"]>=-1e-12
        assert abs(r["balance_residual"])<1e-11


def test_random_both_heat_signs_and_endpoint_Pauli_saturation():
    rng=random.Random(22611)
    signs=set()
    for _ in range(130):
        Ei=rng.uniform(0.,1.);Ej=Ei+rng.uniform(.01,2.)
        fi=rng.uniform(0.,1.);fj=rng.uniform(0.,1.)
        Vi=10**rng.uniform(-1,1);Vj=10**rng.uniform(-1,1)
        A=10**rng.uniform(-1,1);T=10**rng.uniform(-1,1)
        gamma=10**rng.uniform(-1,2);dt=10**rng.uniform(-3,3)
        r=evolve_pair(Ei,Ej,fi,fj,Vi,Vj,A*T**4,A,gamma,dt)
        assert 0<=r["fi"]<=1 and 0<=r["fj"]<=1 and r["ESM"]>=0
        assert abs(r["particle_residual"])<1e-10
        assert abs(r["energy_residual"])<1e-9
        assert r["entropy_delta"]>=-1e-9
        signs.add(1 if r["dN"]>0 else -1)
    assert signs=={-1,1}


def test_network_split_conserves_even_when_time_error_nonzero():
    E=[.1,.6,1.3];V=[1,3,5];f=[.9,.1,.1]
    edges=[(0,1,2.),(1,2,1.)]
    for steps in (1,10,100):
        r=evolve_network(E,V,f,50.,10.,edges,1.,steps)
        assert all(0<=x<=1 for x in r["f"])
        assert abs(r["number_residual"])<1e-12
        assert abs(r["total_energy_residual"])<1e-12
        assert r["entropy_delta"]>=0
    coarse=evolve_network(E,V,f,50.,10.,edges,1.,10)
    fine=evolve_network(E,V,f,50.,10.,edges,1.,100)
    assert max(abs(a-b) for a,b in zip(coarse["f"],fine["f"]))>.001


def test_exact_occupancy_endpoints_return_interior_without_clipping():
    for fi,fj in [(1.,0.),(0.,1.)]:
        r=evolve_pair(.1,1.,fi,fj,1.,1.,10.,1.,1.,100.)
        assert 0<r["fi"]<1 and 0<r["fj"]<1
        assert r["entropy_delta"]>0


def test_zero_portal_and_bad_domain_fail_closed():
    assert evolve_pair(.1,1.,.4,.5,2.,2.,10.,2.,0.,1.)["dN"]==0
    with pytest.raises(ValueError):evolve_pair(.1,.1,.4,.5,2.,2.,10.,2.,1.,1.)
    with pytest.raises(ValueError):evolve_pair(.1,1.,1.1,.5,2.,2.,10.,2.,1.,1.)
    with pytest.raises(ValueError):evolve_pair(.1,1.,.4,.5,2.,2.,10.,2.,-1.,1.)
    with pytest.raises(ValueError):evolve_network([.1,.6],[1,1],[.9,.1],10.,2.,[(1,0,2.)],1.)
