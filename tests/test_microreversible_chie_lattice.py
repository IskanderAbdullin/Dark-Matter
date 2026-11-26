"""Independent physical q-integral χe rate, full Pauli and time-reversible kinetics."""
import math,sys
from pathlib import Path
import numpy as np
import pytest
from scipy.special import expit
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from microreversible_chie_lattice import (
  invariant_kernel,microreversible_lattice,sources,evolve,invariants)


def test_exact_time_reversal_of_independently_evaluated_invariant_rates():
    M=.1;G=.001
    for E,k,Ep,kp in [(1.,2.,2.,1.),(3.,4.,4.,3.),
                       (2.,4.,1.,5.),(.7,1.2,1.4,.5)]:
        fw=float(invariant_kernel(M,G,E,k,np.array([Ep]),np.array([kp]))[0])
        bw=float(invariant_kernel(M,G,Ep,kp,np.array([E]),np.array([k]))[0])
        assert fw>0 and math.isclose(fw,bw,rel_tol=1e-12)
    with pytest.raises(ValueError):
        invariant_kernel(.1,.001,1.,1.,np.array([1.5]),np.array([1.]))


def test_invariant_q_quadrature_and_real_scattering_kernel():
    args=(.1,.001,1.,2.,np.array([2.]),np.array([1.]))
    a=float(invariant_kernel(*args,nq=8)[0])
    b=float(invariant_kernel(*args,nq=16)[0])
    assert abs(a/b-1)<5e-7
    assert a>0
    z=invariant_kernel(.1,0.,1.,2.,np.array([2.]),np.array([1.]))
    assert np.all(z==0)


def test_exact_discrete_gauge_species_number_and_energy_conservation():
    g=microreversible_lattice(nchi=20,ne=22,nq=10)
    E=g["Echi"];k=g["Ee"]
    np.testing.assert_allclose(E[g["i"]]+k[g["a"]],
                               E[g["j"]]+k[g["b"]],rtol=1e-12)
    assert np.all(g["W"]>0)
    assert math.isclose(sum(g["Vchi"]),4/(6*math.pi**2)*
        (((.1+20*.5)**2-.1**2)**1.5),rel_tol=1e-14)
    assert math.isclose(sum(g["Ve"]),4/(6*math.pi**2)*(22*.5)**3,rel_tol=1e-14)
    r=sources(g,expit(-E/.9),expit(-k))
    assert r["Qchi_GeV5"]>0
    assert abs(r["Nchi_source_GeV4"])<1e-17
    assert abs(r["Ne_source_GeV4"])<1e-17
    assert abs(r["energy_residual_GeV5"])<1e-17
    assert r["entropy_GeV4"]>0


def test_exact_common_T_stationary_at_independent_mu_for_entire_source():
    g=microreversible_lattice(nchi=20,ne=20)
    for muc,mue in [(0.,0.),(.3,-.2),(-.1,.2)]:
        E=g["Echi"];k=g["Ee"]
        r=sources(g,expit(muc-E),expit(mue-k))
        assert max(abs(r["dchi"]))<1e-18
        assert max(abs(r["de"]))<1e-18
        assert abs(r["Qchi_GeV5"])<1e-17
        assert abs(r["entropy_GeV4"])<1e-24
        assert r["no_posthoc_rate_symmetrisation"]


def test_independent_unbinned_five_dimensional_Q_and_dE_refinement():
    # Previous independent numerical Gauss 5D weak source, p/k to 22GeV,
    # n28 on p,k,incoming angle,CM-polar and azimuth, distinct kinematics.
    cold_ref=1.9851492233308662e-8
    hot_ref=-4.432408224385169e-8
    for h,n in [(.5,44),(.25,88)]:
        g=microreversible_lattice(dE=h,nchi=n,ne=n,nq=8)
        E=g["Echi"];k=g["Ee"]
        cold=sources(g,expit(-E/.9),expit(-k))["Qchi_GeV5"]
        hot=sources(g,expit(-E/1.1),expit(-k))["Qchi_GeV5"]
        assert abs(cold/cold_ref-1)<.002
        assert abs(hot/hot_ref-1)<.002
        if h==.5:
            assert abs(cold/cold_ref-1)<.0005
            assert abs(hot/hot_ref-1)<.0005


def test_bounded_all_quartet_time_dynamics_with_original_microreversible_W():
    g=microreversible_lattice(dE=.5,nchi=14,ne=17,nq=10)
    E=g["Echi"];k=g["Ee"]
    f=expit(-E/1.2);e=expit(-k)
    Q=sources(g,f,e)["Qchi_GeV5"]
    r=evolve(g,f,e,1e3)
    assert r["Qchi_integrated_GeV4"]<0
    assert abs(r["Qchi_integrated_GeV4"]/(1e3*Q)-1)<.003
    assert abs(r["Nchi_residual"])<1e-12
    assert abs(r["Ne_residual"])<1e-12
    assert abs(r["energy_residual_GeV4"])<1e-12
    assert r["entropy_gain_GeV3"]>0
    assert np.all((r["fchi"]>=0)&(r["fchi"]<=1))
    assert np.all((r["fe"]>=0)&(r["fe"]<=1))
    assert r["edge_count"]>900


def test_stationary_FD_under_stiff_all_quartet_time():
    g=microreversible_lattice(dE=.5,nchi=14,ne=17,nq=10)
    E=g["Echi"];k=g["Ee"]
    f=expit(.15-E);e=expit(-.2-k)
    r=evolve(g,f,e,1e6,substeps=4)
    assert max(abs(r["fchi"]-f))<1e-11
    assert max(abs(r["fe"]-e))<1e-11
    assert abs(r["entropy_gain_GeV3"])<1e-10


@pytest.mark.parametrize("kwargs",[{"G":-1.},{"M":0.},{"nchi":2},
                                    {"ne":2},{"dE":0.},{"nq":1}])
def test_bad_inputs_fail_closed(kwargs):
    with pytest.raises(ValueError):
        microreversible_lattice(**kwargs)
