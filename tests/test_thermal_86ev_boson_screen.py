"""Independent conditional hot/thermal real 86-eV boson cosmology screen.

These tests deliberately do not claim the full Yukawa thermal plasma produces
this spectrum, nor translate a fermionic Lyman-alpha mass bound into a boson
95% exclusion without transfer-function/likelihood calculations.
"""
import math

from cosmology.thermal_86ev_boson_screen import (
    thermal_abundance,thermal_momentum_weights,
    temperature_ratio_for_total_dark_matter,mean_collisionless_horizon,
    T_GAMMA_K,EV_PER_K,OMEGA_C_H2
)


def test_PDG_entropy_ratio_and_bosonic_hot_relic_abundance():
    state=thermal_abundance(86.,61.75)
    assert math.isclose(state["Tphi_over_Tgamma"],
                        .3985439640892993,rel_tol=3e-12)
    assert math.isclose(state["nphi_cm_minus3"],
                        12.999403238866398,rel_tol=3e-12)
    assert math.isclose(state["Omega_phi_h2"],
                        .10610025496952659,rel_tol=3e-12)
    assert math.isclose(state["delta_Neff_early_thermal"],
                        .05554485291636691,rel_tol=3e-12)
    assert math.isclose(state["gstar_entropy_dec_for_omega_c_h2"],
                        54.597422869735546,rel_tol=3e-12)
    for g in (15.,61.75,106.75):
        s=thermal_abundance(86.,g)
        # For the free adiabatic thermal relic, Ωh²∝1/g*S,dec.
        assert math.isclose(s["Omega_phi_h2"]*g,
                            86.*410.7*.5*3.909/10536.72,rel_tol=3e-12)


def test_independent_exact_number_weighted_Bose_Fermi_momentum_distributions():
    for kind,expected_x in (("real_boson",2.7011780329190644),
                            ("two_fermion",3.1513743717389078)):
        x,w=thermal_momentum_weights(kind,96)
        assert len(x)==len(w)==96
        assert all(v>0 for v in w)
        assert math.isclose(sum(w),1.,rel_tol=2e-13)
        assert math.isclose(sum(w*x),expected_x,rel_tol=3e-12)


def test_86eV_all_DM_scalar_has_order_Mpc_characteristic_free_streaming():
    r=temperature_ratio_for_total_dark_matter(86.,"real_boson")
    assert math.isclose(r,.4152387187558684,rel_tol=3e-12)
    t=r*T_GAMMA_K*EV_PER_K
    results=[mean_collisionless_horizon(86.,t,"real_boson",n)
             for n in (48,96,128)]
    for a in results:
        assert .0103<a["mean_velocity_eq_over_c"]<.0106
        assert 6.5<a["mean_comoving_freestream_to_eq_Mpc"]<6.8
    assert abs(results[1]["mean_comoving_freestream_to_eq_Mpc"]
               -results[2]["mean_comoving_freestream_to_eq_Mpc"])<2e-5
    assert math.isclose(results[-1]["mean_comoving_freestream_to_eq_Mpc"],
                        6.6400167562020656,rel_tol=3e-6)


def test_same_abundance_thermal_FD_5p7keV_is_much_colder_but_not_same_lyman_likelihood():
    ratio=temperature_ratio_for_total_dark_matter(5700.,"two_fermion")
    assert .0895<ratio<.0898
    w=mean_collisionless_horizon(5700.,ratio*T_GAMMA_K*EV_PER_K,
                                 "two_fermion",128)
    assert math.isclose(w["mean_comoving_freestream_to_eq_Mpc"],
                        .05547491887238949,rel_tol=3e-6)
    assert 3.9e-5<w["mean_velocity_eq_over_c"]<4.1e-5
    r=temperature_ratio_for_total_dark_matter(86.,"real_boson")
    b=mean_collisionless_horizon(86.,r*T_GAMMA_K*EV_PER_K,
                                 "real_boson",128)
    assert b["mean_comoving_freestream_to_eq_Mpc"]>100*(
        w["mean_comoving_freestream_to_eq_Mpc"])
    # Pure WDM Lyman-alpha 5.7-keV is external observational context;
    # this path-length ratio is NOT a transferred 95% C.L. constraint.


def test_stepfunction_gstar_and_collisionless_horizon_are_assumptions_not_yukawa_prediction():
    early=thermal_abundance(86.,61.75)
    allDM=thermal_abundance(86.,early["gstar_entropy_dec_for_omega_c_h2"])
    assert math.isclose(allDM["Omega_phi_h2"],OMEGA_C_H2,rel_tol=3e-13)
    r=mean_collisionless_horizon(86.,early["Tphi0_eV"],"real_boson",128)
    assert math.isclose(r["mean_comoving_freestream_to_eq_Mpc"],
                        6.428442647590668,rel_tol=3e-6)
    assert r["mean_velocity_eq_over_c"]>.009


def test_analytic_path_crosscheck_existing_verified_CLASS_86eV_real_boson_gstar55():
    # CLASS was already run on 2026-11-21, see THERMAL-CLASS-VERIFIED.md.
    a=thermal_abundance(86.,55.)
    assert math.isclose(a["Omega_phi_h2"],.11912164989760488,rel_tol=3e-12)
    # The pre-existing CLASS analytic baseline used a more precise photon
    # density than the present deliberately rounded n_gamma=410.7/cm^3;
    # do not require bitwise-identical Ω from distinct stated inputs.
    assert abs(a["Omega_phi_h2"]/.11912939486124709-1)<7e-5
    q=mean_collisionless_horizon(86.,a["Tphi0_eV"],"real_boson",128)
    assert math.isclose(q["mean_comoving_freestream_to_eq_Mpc"],6.627211308933845,rel_tol=3e-6)
    khalf_existing_CLASS=.43275181165
    assert math.isclose(khalf_existing_CLASS*q["mean_comoving_freestream_to_eq_Mpc"],2.867937700128489,rel_tol=3e-6)
    # This product is a consistency diagnostic, NOT a universal mapping.
