"""All-edge two-species 5D-derived Pauli lattice dynamics and order error."""
import sys
from pathlib import Path
import numpy as np
import pytest
from scipy.special import expit
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from paired_energy_quartets import quartets,sources
from paired_quartet_network import evolve,edges_from_paired_W


def model():
    return quartets(nchi=14,ne=17,dE=.5,nmu=8,ncm=8,nphi=8)


def test_all_quartets_dynamically_evolve_chi_and_e_without_reservoir():
    g=model()
    assert len(edges_from_paired_W(g))>1000
    f=expit(-g["Echi"]/1.2);e=expit(-g["Ee"])
    for steps in (1,4,16):
        r=evolve(g,f,e,2e6,substeps=steps)
        assert r["delta_Echi"]<0 and r["delta_Ee"]>0
        assert abs(r["Nchi_residual"])<1e-12
        assert abs(r["Ne_residual"])<1e-12
        assert abs(r["total_energy_residual"])<1e-12
        assert r["entropy_gain"]>0
        assert np.all((r["fchi"]>=0)&(r["fchi"]<=1))
        assert np.all((r["fe"]>=0)&(r["fe"]<=1))
        assert not r["is_original_unmodified_5D_kernel"]


def test_common_T_independent_mu_is_stationary_even_at_large_dt():
    g=model();f=expit(.15-g["Echi"]);e=expit(-.2-g["Ee"])
    r=evolve(g,f,e,1e7,substeps=4)
    assert max(abs(r["fchi"]-f))<1e-12
    assert max(abs(r["fe"]-e))<1e-12
    assert abs(r["entropy_gain"])<1e-12


def test_short_time_slope_agrees_with_independent_instantaneous_source():
    g=model();f=expit(-g["Echi"]/1.2);e=expit(-g["Ee"])
    Q=sources(g,f,e)["Qchi"]
    for dt,limit in [(1e3,.001),(1e4,.005)]:
        r=evolve(g,f,e,dt,substeps=1)
        assert abs(r["delta_Echi"]/dt/Q-1)<limit


def test_commutator_order_error_is_real_and_decreases_on_dt_refinement():
    g=model();f=expit(-g["Echi"]/1.2);e=expit(-g["Ee"])
    deviations=[]
    for steps in (1,4,16):
        a=evolve(g,f,e,2e6,substeps=steps)
        b=evolve(g,f,e,2e6,substeps=steps,reverse=True)
        assert abs(a["total_energy_residual"])<1e-12
        assert abs(b["total_energy_residual"])<1e-12
        deviations.append(abs(a["delta_Echi"]-b["delta_Echi"]))
    assert deviations[0]>1e-4
    assert deviations[2]<deviations[1]<deviations[0]
    assert deviations[2]<1e-5


def test_reject_bad_time_or_occupancies():
    g=model();f=expit(-g["Echi"]);e=expit(-g["Ee"])
    with pytest.raises(ValueError):evolve(g,f,e,-1)
    with pytest.raises(ValueError):evolve(g,f,e,1,substeps=0)
    with pytest.raises(ValueError):evolve(g,-f,e,1)
    with pytest.raises(ValueError):evolve(g,f,e[:-1],1)
