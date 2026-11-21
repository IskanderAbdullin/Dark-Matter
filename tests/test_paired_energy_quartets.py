"""Numerical microreversible Pauli χ+e quartet checks and negative controls."""
import sys
import math
from pathlib import Path
import numpy as np
import pytest
from scipy.special import expit
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from paired_energy_quartets import quartets,sources,implicit_one_quartet


def test_physical_capacities_and_exact_quartet_energy():
    x=quartets(nchi=18,ne=26,dE=.5,nmu=8,ncm=8,nphi=8)
    assert math.isclose(sum(x["Vchi"]),4/(6*math.pi**2)*
        (((18*.5+.1)**2-.1**2)**1.5),rel_tol=1e-14)
    assert math.isclose(sum(x["Ve"]),4/(6*math.pi**2)*(26*.5)**3,rel_tol=1e-14)
    i,a,j,b=np.nonzero(x["W"]>0)
    np.testing.assert_allclose(x["Echi"][i]+x["Ee"][a],
        x["Echi"][j]+x["Ee"][b],rtol=1e-13,atol=1e-12)
    assert np.all(x["W"]>=0)
    assert np.array_equal(x["W"],x["W"].transpose(2,3,0,1))


def test_common_T_FD_stops_for_independent_chemical_potentials():
    x=quartets(nchi=18,ne=26,dE=.5,nmu=8,ncm=8,nphi=8)
    E=x["Echi"];K=x["Ee"]
    for muc,mue in [(0.,0.),(.2,-.1),(-.2,.3)]:
        r=sources(x,expit(muc-E),expit(mue-K))
        assert abs(r["Qchi"])<1e-17 and abs(r["Qe"])<1e-17
        assert max(abs(r["dchi"]))<1e-17
        assert max(abs(r["de"]))<1e-17
        assert abs(r["entropy"])<1e-25


def test_hotter_chi_cools_colder_chi_heats_and_total_Q_zero():
    x=quartets(nchi=36,ne=54,dE=.5,nmu=10,ncm=10,nphi=10)
    E=x["Echi"];K=x["Ee"];electron=expit(-K)
    for Tc,sign in [(.9,+1),(1.1,-1)]:
        r=sources(x,expit(-E/Tc),electron)
        assert sign*r["Qchi"]>0
        assert abs(r["Qe"]+r["Qchi"])<1e-18
        assert abs(r["number_chi"])<1e-18
        assert abs(r["number_e"])<1e-18
        assert r["entropy"]>0
    assert x["W_eq_Pauli_asym_L1"]>.01
    assert abs(sources(x,raw=True)["Qchi"])>1e-13


def test_independent_full_5D_unbinned_weak_heat_comparison():
    x=quartets(nchi=36,ne=54,dE=.5,nmu=10,ncm=10,nphi=10)
    E=x["Echi"];K=x["Ee"];fe=expit(-K)
    # Different continuous 5D p/k Gauss implementation, n14 in each
    # dimension, momenta p,k in [0,22GeV], no discrete quartet projection.
    cold_ref=1.9852809333031555e-8
    hot_ref=-4.432526721786192e-8
    assert abs(sources(x,expit(-E/.9),fe)["Qchi"]/cold_ref-1)<.025
    assert abs(sources(x,expit(-E/1.1),fe)["Qchi"]/hot_ref-1)<.025
    assert x["omitted_FD_weight_fraction"]<1e-4


def test_physical_rate_single_quartet_stiff_implicit_Pauli_H():
    x=quartets(nchi=30,ne=30,dE=.4,nmu=8,ncm=8,nphi=8)
    E=x["Echi"];K=x["Ee"];W=x["W"];f=expit(-E);g=expit(-K)
    I,A,J,B=np.unravel_index(
        np.argmax(W*(E[:,None,None,None]<E[None,None,:,None])),W.shape)
    for dt in (100.,1e5,1e8):
        r=implicit_one_quartet(x,int(I),int(A),int(J),int(B),
            float(.7*f[I]),float(1.1*g[A]),float(1.1*f[J]),
            float(.9*g[B]),dt)
        assert all(0<=r[t]<=1 for t in ("fi","ga","fj","gb"))
        assert abs(r["chi_charge_residual"])<1e-10
        assert abs(r["electron_charge_residual"])<1e-10
        assert abs(r["energy_residual_GeV4"])<1e-12
        assert r["entropy_change_GeV3"]>=0
        assert r["not_exact_unprojected_microphysics"]


def test_raw_asymmetry_reduces_with_angular_refinement():
    a=quartets(nchi=24,ne=24,dE=.5,nmu=8,ncm=8,nphi=8)
    b=quartets(nchi=24,ne=24,dE=.5,nmu=16,ncm=16,nphi=16)
    assert b["W_eq_Pauli_asym_L1"]<a["W_eq_Pauli_asym_L1"]
    assert a["W_eq_Pauli_asym_L1"]>.05


@pytest.mark.parametrize("args",[{"M":0.},{"G":-1.},{"nchi":2},
                                  {"dE":0.},{"ne":4000},{"ge":0.}])
def test_fail_closed(args):
    with pytest.raises(ValueError):quartets(**args)
