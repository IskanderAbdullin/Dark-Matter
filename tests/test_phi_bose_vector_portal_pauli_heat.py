"""Finite vector A' BOTH Dirac χ1/χ2 and all antiparticles, physical
Pauli 5D electron weak heat; independent contact 5D reference and entropy.

All electron/e+ thermal bath results below are CONDITIONAL partial SM:
no photon/nu/hadrons, finite electron mass, in-medium mediator, strong
coupled SM, numerical drag gamma_el, FRW or 86-eV relic.
"""
import math

from kinetics.relativistic_weak_form_2to2 import weak_energy_moments
from kinetics.phi_bose_vector_portal_finite_mediator import vector_portal_setup
from kinetics.phi_bose_vector_portal_pauli_heat import (
    one_dirac_charge_sector_weak_heat,four_dirac_sectors_electron_heat,
)


def test_independent_physical_finite_mediator_weak_heat_reduces_to_existing_contact_5D():
    # The old module constructs its CM boost from the exact lab frame
    # and separately integrates the Pauli forward-reverse weak energy.
    m=.02
    # Genuine asymptotic contact check: mV=1000 GeV still has a finite
    # 6.87e-5 high-energy-tail propagator correction on this cutoff grid;
    # use 10 TeV to make the contact error <2e-6 instead of hiding it.
    massive=10000.
    coupling=vector_portal_setup(mV=massive)["G_contact_GeV_minus2"]
    for n in (6,8,10):
        finite=one_dirac_charge_sector_weak_heat(
            mchi=m,Tchi=.8,Te=1.,mV=massive,
            n_p=n,n_k=n,n_mu=n,n_cm=n,n_az=n,cutoff=18.)
        contact=weak_energy_moments(
            M=m,Tchi=.8,Te=1.,G=coupling,
            gchi=2,ge=4,n_p=n,n_k=n,n_mu=n,n_cm=n,n_phi=n,
            cutoff=18.)
        assert finite["Q_chi_GeV5"]>0
        assert contact["Q_chi_GeV5"]>0
        assert math.isclose(finite["Q_chi_GeV5"],
                            contact["Q_chi_GeV5"],rel_tol=3e-5)
        assert math.isclose(finite["quantum_H_GeV4"],
                            contact["entropy_GeV4"],rel_tol=3e-5)
        assert abs(finite["Q_chi_GeV5"]+finite["Q_e_GeV5"])<1e-30


def test_independent_five_dimensional_weak_heat_positive_and_reverse_sign():
    for m in (.02,.1):
        for temp in (.8,1.,1.2):
            v=one_dirac_charge_sector_weak_heat(
                mchi=m,Tchi=temp,Te=1.,n_p=8,n_k=8,
                n_mu=8,n_cm=8,n_az=8,cutoff=18.)
            assert v["max_four_shell_relative_defect"]<1e-12
            assert v["microscopic_nonnegative_H"]
            assert v["quantum_H_GeV4"]>=0
            assert v["total_positive_event_measure"]>0
            assert abs(v["Q_chi_GeV5"]+v["Q_e_GeV5"])<1e-30
            if temp<1:
                assert v["Q_chi_GeV5"]>0 and v["quantum_H_GeV4"]>0
            if temp>1:
                assert v["Q_chi_GeV5"]<0 and v["quantum_H_GeV4"]>0
            if temp==1.:
                assert v["Q_chi_GeV5"]==0 and v["quantum_H_GeV4"]==0


def test_full_both_dirac_and_anti_dark_electron_heat_exact_conservation_and_H():
    r=four_dirac_sectors_electron_heat(Tc1=.8,Tc2=.9,Te=1.,n=8)
    assert len(r["components"])==4
    assert {x["species"] for x in r["components"]}=={
        "chi1","anti_chi1","chi2","anti_chi2"}
    assert r["Q_dark_from_e_GeV5"]>0
    assert r["Q_e_bath_GeV5"]<0
    assert abs(r["exact_heat_ledger_Qdark_plus_Qe_GeV5"])<1e-28
    assert r["H_total_GeV4"]>0
    assert all(x["Q_chi_GeV5"]>0 for x in r["components"])
    for i,j in ((0,1),(2,3)):
        assert math.isclose(r["components"][i]["Q_chi_GeV5"],
                            r["components"][j]["Q_chi_GeV5"],rel_tol=1e-13)


def test_physical_finite_vector_propagator_is_not_equivalent_contact_at_1GeV_heat():
    m=.02
    for n in (6,8,10):
        v=one_dirac_charge_sector_weak_heat(
            mchi=m,Tchi=.8,Te=1.,mV=10.,
            n_p=n,n_k=n,n_mu=n,n_cm=n,n_az=n,cutoff=18.)
        G=vector_portal_setup(mV=10.)["G_contact_GeV_minus2"]
        c=weak_energy_moments(M=m,Tchi=.8,Te=1.,G=G,
                              gchi=2,ge=4,n_p=n,n_k=n,
                              n_mu=n,n_cm=n,n_phi=n,cutoff=18.)
        assert 0<v["Q_chi_GeV5"]<c["Q_chi_GeV5"]
        ratio=v["Q_chi_GeV5"]/c["Q_chi_GeV5"]
        assert .35<ratio<1.
        print("VECTOR_FINITE_MASS_PAULI_HEAT",repr(dict(
            n=n,Q_fin=v["Q_chi_GeV5"],Q_contact=c["Q_chi_GeV5"],
            finite_to_contact_ratio=ratio,H=v["quantum_H_GeV4"])),flush=True)


def test_vector_portal_weak_heat_scales_couplings_squared_and_charged_species_degree():
    strong=one_dirac_charge_sector_weak_heat(
        mchi=.02,Tchi=.8,Te=1.,n_p=6,n_k=6,
        n_mu=6,n_cm=6,n_az=6,cutoff=18.)
    weak=one_dirac_charge_sector_weak_heat(
        mchi=.02,Tchi=.8,Te=1.,epsilon=5e-6,
        n_p=6,n_k=6,n_mu=6,n_cm=6,n_az=6,cutoff=18.)
    assert math.isclose(weak["Q_chi_GeV5"]/strong["Q_chi_GeV5"],.25,rel_tol=1e-12)
    assert math.isclose(weak["quantum_H_GeV4"]/strong["quantum_H_GeV4"],.25,rel_tol=1e-12)


def test_full_vector_weak_chemical_potential_is_not_sm_bath_chemical_equilibration():
    a=one_dirac_charge_sector_weak_heat(
        mchi=.1,Tchi=1.,Te=1.,mu_chi=.03,mu_e=-.02,
        n_p=6,n_k=6,n_mu=6,n_cm=6,n_az=6,cutoff=18.)
    assert a["Q_chi_GeV5"]==0
    assert a["quantum_H_GeV4"]==0
    assert a["number_chi_source_GeV4"]==0
    assert a["number_e_source_GeV4"]==0
    # Purely elastic vector events preserve each number and are stationary
    # at common T even when the two FD chemical potentials differ.
