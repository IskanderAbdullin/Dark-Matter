"""Detailed balance and dimensional matching, not external BEST execution."""
import math
import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from microphysics_thermal_gain_loss import (
    fermi,bose,gain_loss_on_shell,yukawa_best_bridge)
from microphysics_portal import width_from_best_spin_sum


@pytest.mark.parametrize("T,mu",[(1.,0.),(.3,.15),(5.,-.2)])
def test_fermi_bose_ratio_identities(T,mu):
    E=1.6
    f=fermi(E,T,mu)
    b=bose(E,T,mu)
    exp=math.exp(-(E-mu)/T)
    assert math.isclose(f/(1-f),exp,rel_tol=1e-13)
    assert math.isclose(b/(1+b),exp,rel_tol=1e-13)


@pytest.mark.parametrize("Ephi,T",[(.4,.3),(.01,.2),(1.5,1.)])
def test_chemical_equilibrium_with_real_boson_zero_chemical_potential(Ephi,T):
    q=gain_loss_on_shell(1.1+Ephi,1.1,Ephi,T,mu1=.1,mu2=.1,muphi=0.)
    assert abs(q["net_factor"])<3e-14*max(1.,q["decay_gain_factor"])
    assert abs(q["detailed_balance_stable_net"])<1e-14
    assert q["is_full_finite_temperature_rate"] is False


def test_signed_affinity_matches_exact_thermal_ratio():
    p=gain_loss_on_shell(1.6,1.1,.5,1.,mu1=.3,mu2=.1)
    n=gain_loss_on_shell(1.6,1.1,.5,1.,mu1=-.1,mu2=.1)
    assert p["net_factor"]>0 and n["net_factor"]<0
    for r in [p,n]:
        assert math.isclose(r["decay_gain_factor"]/r["inverse_gain_factor"],
                            math.exp(r["chemical_affinity_over_T"]),rel_tol=1e-13)
        assert math.isclose(r["net_factor"],r["detailed_balance_stable_net"],
                            rel_tol=1e-13,abs_tol=1e-14)


def test_large_occupation_not_same_as_density_or_BEC():
    a=gain_loss_on_shell(1.001,1.,.001,1.,mu1=0.,mu2=0.)
    assert a["occupations"]["bphi"]>900
    assert abs(a["net_factor"])<1e-12
    assert a["is_full_finite_temperature_rate"] is False


def test_exact_dimensionless_BEST_normalization_bridge():
    z=yukawa_best_bridge(.1,.002,86e-9,.3,2e-8)
    width_recovered=width_from_best_spin_sum(
        z["BEST_m1_over_T"],z["BEST_m2_over_T"],z["BEST_mphi_over_T"],
        z["BEST_spin_summed_m2_over_T2"],4)
    assert math.isclose(width_recovered,z["BEST_Gamma_over_T"],rel_tol=1e-13)
    assert math.isclose(z["physical_Gamma_rest_GeV"],
                        z["BEST_Gamma_over_T"]*.3,rel_tol=1e-13)
    assert not z["external_BEST_executed"]
    assert not z["CLASS_ready"]


def test_thermal_invalid_onshell_and_chemical_potential_fail_closed():
    with pytest.raises(ValueError): gain_loss_on_shell(1.,.1,.2,1.)
    with pytest.raises(ValueError): gain_loss_on_shell(1.,.5,.5,1.,muphi=.5)
    with pytest.raises(ValueError): bose(.1,1.,mu=.11)
    with pytest.raises(ValueError): yukawa_best_bridge(.1,.2,.01,1.,1e-4)
