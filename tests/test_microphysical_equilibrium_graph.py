"""Tests: real 5D Pauli equilibrium rates, bias and fake reversal controls."""
import math
import sys
from pathlib import Path
import numpy as np
import pytest
from scipy.special import expit
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from microphysical_equilibrium_graph import (
    graph,calibrated_pair_rates,graph_instantaneous_heat,
    local_pair_implicit_step,cm_event)
from relativistic_weak_form_2to2 import cm_event as prior_cm


def test_independent_on_shell_cm_kinematic_agreement():
    u=np.array([-.8,-.3,.6])[:,None]
    phi=np.array([0.,1.,3.])[None,:]
    for M,p,k,mu in [(.1,0.,2.,-.3),(.1,.7,3.,.5),(1.,.2,.05,-.8)]:
        ours=cm_event(M,p,k,mu,u,phi)
        ref=prior_cm(M,p,k,mu,u,phi)
        for a,b in zip(ours,ref[:3]):
            np.testing.assert_allclose(a,b,atol=1e-14,rtol=1e-13)
        assert np.max(ref[3])<1e-12


def test_exact_physical_fermion_capacities_and_closed_FD():
    x=graph(nchi=24,nk=12,nmu=10,ncm=10,nphi=10)
    E=x["E"];f=expit(-E)
    assert math.isclose(sum(x["capacities"]),4*22**3/(6*math.pi**2),rel_tol=1e-14)
    assert np.all(x["capacities"]>0)
    assert np.all(x["directed"]>=0)
    assert np.array_equal(x["symmetric"],x["symmetric"].T)
    assert x["raw_stationary_relative_L1"]>.01  # failure BEFORE closure
    assert x["outgoing_out_of_support_weight_fraction"]<2e-6
    Q,net=graph_instantaneous_heat(x,f)
    assert abs(Q)<1e-17 and np.max(np.abs(net))<1e-18
    assert all(rate>=0 for _,_,rate in calibrated_pair_rates(x))


def test_jensen_bias_positive_even_when_event_energy_is_preserved():
    xx=[graph(nchi=n,nk=12,nmu=10,ncm=10,nphi=10) for n in (16,32,64)]
    for x in xx:
        assert x["D_graph"]>x["D_true_inside"]>0
        assert abs(x["D_graph"]-x["D_graph_symmetric"])<1e-20
        assert x["D_true_all"]>=x["D_true_inside"]
    assert xx[0]["jensen_relative"]>.05
    assert xx[2]["jensen_relative"]<.005
    assert xx[0]["jensen_relative"]>xx[1]["jensen_relative"]>xx[2]["jensen_relative"]


def test_continuum_microphysical_5D_energy_susceptibility_and_nonlinear_slope():
    x=graph(nchi=64,nk=16,nmu=16,ncm=16,nphi=16)
    # Independently measured D11 in original 5D FD Dirichlet n24/cut32.
    independent_D11=3.027899802221842e-7
    assert abs(x["D_true_all"]/independent_D11-1)<.003
    assert abs(x["D_graph"]/independent_D11-1)<.008
    assert .01<x["unpaired_asymmetry"]<.03
    E=x["E"]; eps=1e-4
    qp,_=graph_instantaneous_heat(x,expit(-E+eps*E))
    qm,_=graph_instantaneous_heat(x,expit(-E-eps*E))
    assert math.isclose((qp-qm)/(2*eps),-x["D_graph"],rel_tol=1e-6)
    assert qp<0 and qm>0


def test_rejected_unpaired_flux_decreases_with_angular_refinement():
    a=graph(nchi=40,nk=12,nmu=12,ncm=12,nphi=12)
    b=graph(nchi=40,nk=20,nmu=20,ncm=20,nphi=20)
    assert b["unpaired_asymmetry"]<a["unpaired_asymmetry"]
    assert b["raw_stationary_relative_L1"]<a["raw_stationary_relative_L1"]
    assert abs(b["D_true_all"]/a["D_true_all"]-1)<.002


def test_finite_bath_is_conditional_frozen_T0_graph_not_physical_all_T():
    g=graph(nchi=20,nk=10,nmu=8,ncm=8,nphi=8)
    i,j,_=max(calibrated_pair_rates(g),key=lambda r:r[2])
    E=g["E"];f=expit(-E)
    r=local_pair_implicit_step(g,i,j,float(f[i]),float(f[j]),10.,10.,1e6)
    assert abs(r["dN"])<1e-9
    fi=float(expit(-E[i]+.02*E[i]));fj=float(expit(-E[j]+.02*E[j]))
    v=local_pair_implicit_step(g,i,j,fi,fj,10.,10.,1e6)
    assert abs(v["particle_residual"])<1e-12
    assert abs(v["energy_residual"])<1e-12
    assert v["entropy_delta"]>=0
    assert 0<=v["fi"]<=1 and 0<=v["fj"]<=1
    assert v["physical_nonthermal_full_collision"] is False
    assert abs(v["T"]-1)<.001


@pytest.mark.parametrize("kwargs",[{"M":0.},{"G":-1.},{"nchi":2},
                                    {"nk":2},{"pmax":0.}])
def test_reject_nonphysical_inputs(kwargs):
    with pytest.raises(ValueError):graph(**kwargs)


def test_zero_portal_is_well_defined_without_spurious_division():
    x=graph(G=0.,nchi=6,nk=4,nmu=3,ncm=4,nphi=4)
    assert x["D_graph"]==x["D_true_all"]==0.
    assert x["unpaired_asymmetry"]==0.
    assert x["raw_FD_heat_without_reversible_closure_GeV5"]==0.
    assert calibrated_pair_rates(x)==[]


def test_raw_microprojected_FD_heat_is_false_before_reversible_closure():
    x=graph(nchi=24,nk=12,nmu=10,ncm=10,nphi=10)
    # At exactly equal χ/e temperatures the microscopic continuum heat is 0.
    # Finite-bin quadrature creates a false signal, removed only by a
    # declared model-changing symmetrization of equilibrium conductances.
    assert x["raw_FD_heat_without_reversible_closure_GeV5"]>1e-9
    assert abs(x["raw_FD_heat_without_reversible_closure_GeV5"]/
               x["D_true_all"])>.003
    q,_=graph_instantaneous_heat(x,expit(-x["E"]))
    assert abs(q)<1e-17
