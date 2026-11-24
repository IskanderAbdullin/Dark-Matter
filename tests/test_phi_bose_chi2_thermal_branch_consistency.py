"""Physical χ2 Yukawa exact thermal pair/H vs tiny electron transport/H."""
import math
from kinetics.phi_bose_chi2_thermal_branch_consistency import (
    instantaneous_thermal_consistency,
)


def test_declared_vector_portal_does_not_ensure_sm_thermal_chi2_at_late_freezeout():
    rows=[instantaneous_thermal_consistency(x) for x in (3,5,8,10,12,15)]
    for r in rows:
        assert r["chi2_to_phi_phi_full_GG_sigma_v_GeV_minus2"]>0
        assert r["chi2_pair_forward_per_chi_over_H"]>0
        assert r["finite_mediator_e_transport_proxy_over_H_NOT_drag"]>0
        assert r["chi2_pair_offshell_vector_e_per_chi_over_H"]>0
        assert r["plasma_electron_assumed_massless"]
        assert r["no_physical_gamma_el_or_full_SM_FWR"]
    print("CHI2_THERMAL_ASSUMPTION_GATE",repr([dict(
            x=r["x"],T_MeV=r["T_GeV"]*1000.,
            pairGamma_over_H=r["chi2_pair_forward_per_chi_over_H"],
            electron_transport_over_H=r[
                "finite_mediator_e_transport_proxy_over_H_NOT_drag"],
            offshell_vector_pair_to_e_over_H=r[
                "chi2_pair_offshell_vector_e_per_chi_over_H"])
                for r in rows]),flush=True)
    assert 2<rows[2]["chi2_pair_forward_per_chi_over_H"]<10
    assert rows[2]["finite_mediator_e_transport_proxy_over_H_NOT_drag"]<1e-6
    assert rows[3]["finite_mediator_e_transport_proxy_over_H_NOT_drag"]<1e-6
    assert rows[4]["chi2_pair_forward_per_chi_over_H"]<.15
    assert all(r["chi2_pair_offshell_vector_e_per_chi_over_H"]<
               1e-7*r["chi2_pair_forward_per_chi_over_H"] for r in rows)


def test_full_GG_thermal_pair_rate_consistent_with_NR_asymptotics_but_not_x3():
    a=instantaneous_thermal_consistency(3.)
    b=instantaneous_thermal_consistency(20.)
    assert math.isclose(a["exact_GG_to_pwave_6b_over_x"],
                        .4813328853691621,rel_tol=2e-8)
    assert math.isclose(b["exact_GG_to_pwave_6b_over_x"],
                        .888605994,rel_tol=5e-7)
    assert a["no_actual_relativistic_relic_prediction"] is True
    assert b["dark_phi_thermal_contact_not_proven"]
