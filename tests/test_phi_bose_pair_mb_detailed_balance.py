"""Vacuum pair detailed balance needs BOTH final 1/2 and inverse initial 1/2.

Cross-checks an independently computed inverse total cross section and
physical Maxwell-Boltzmann thermal forward/inverse reaction integrals.
This is NOT a quantum FD/BE evolving physical FRW source.
"""
import math
import numpy as np
from kinetics.phi_bose_pair_annihilation_mb import (
    M2,MPHI,Y,cross_section,thermal_forward_rate_MB,
)
from kinetics.phi_bose_pair_mb_detailed_balance import (
    reverse_cross_section,equilibrium_pair_densities,
    thermal_MB_chemical_source,
)


def test_spin_and_identical_pair_inverse_cross_section_physical_flux_reciprocity():
    for energy in (.041,.05,.1,.5,1.,5.,100.):
        s=energy**2
        forward=cross_section(s,nangle=96)
        backward=reverse_cross_section(s,nangle=96)
        assert forward>0 and backward>0
        pchi2=s/4-M2*M2
        pphi2=s/4-MPHI*MPHI
        assert math.isclose(backward/forward,
                            8*pchi2/pphi2,rel_tol=3e-14)
        # Exact invariant common 2->2 s event measure:
        # gχ*ganti=4 vs gφ²/(2!)=.5. Sum spin ×4,
        # removal of identical-final 1/2: ×2.
        assert math.isclose(
            .5*backward*(s-4*MPHI*MPHI),
            4*forward*(s-4*M2*M2),rel_tol=5e-14)
    assert reverse_cross_section((2*M2)**2)==0


def test_independent_forward_reverse_MB_reaction_densities_cancel_thermally():
    for T,gstar in ((.1,17.25),(1.,61.75)):
        q=equilibrium_pair_densities(T,nangle=48)
        direct=thermal_forward_rate_MB(T,gstar,nangle=48)
        assert q["gamma_forward_equilibrium_GeV4"]>0
        assert math.isclose(
            q["gamma_forward_equilibrium_GeV4"],
            direct["forward_pair_reaction_density_GeV4"],
            rel_tol=2e-5)  # independent adaptive s-integrals: epsrel 3-5e-6
        assert abs(q["inverse_to_forward_ratio"]-1)<2e-5
        assert math.isclose(
            q["gamma_forward_equilibrium_GeV4"],
            q["gamma_inverse_equilibrium_GeV4"],rel_tol=2e-5)
        assert q["average_total_event_energy_GeV"]>2*M2


def test_common_event_measure_yields_zero_equilibrium_net_flux_and_H():
    for mu2 in (-.01,0.,+.15):
        out=thermal_MB_chemical_source(1.,mu2,-mu2,0.,nangle=48)
        assert out["net_events_GeV4"]==0.
        assert out["H_per_volume_GeV4"]==0.
        assert out["dn_phi_GeV4"]==0.
        assert out["total_energy_defect_GeV5"]==0.
        assert out["net_Dirac_U1_defect"]==0.
        assert abs(sum((out[k] for k in
                        ("dE_chi_GeV5","dE_anti_GeV5","dE_phi_GeV5"))))<1e-30


def test_pair_MB_affinity_counts_energy_and_entropy_positive_in_both_directions():
    for potentials,sign in (((.05,-.05,-.2),+1),
                            ((.05,-.05,+.1),-1),
                            ((+.25,-.05,0.),+1)):
        z=thermal_MB_chemical_source(1.,*potentials,nangle=48)
        assert sign*z["net_events_GeV4"]>0
        assert z["H_per_volume_GeV4"]>0
        assert math.isclose(z["dn_chi_GeV4"],-z["net_events_GeV4"])
        assert math.isclose(z["dn_phi_GeV4"],+2*z["net_events_GeV4"])
        assert z["dn_chi_GeV4"]-z["dn_anti_GeV4"]==0.
        assert math.isclose(z["dE_phi_GeV5"],
                            -z["dE_chi_GeV5"]-z["dE_anti_GeV5"])
        assert math.isclose(z["H_per_volume_GeV4"],
                            z["net_events_GeV4"]*
                            z["chemical_affinity_dimensionless"],rel_tol=2e-13)


def test_physical_pair_total_energy_derived_from_beta_derivative_independently():
    # With a true MB event integral γ(T)=∫W exp[-β(E1+E2)],
    # energy flux ∫W(E1+E2)exp(-βE) = -dγ/dβ = T²dγ/dT.
    T=1.
    h=.004
    low=thermal_forward_rate_MB(T-h,61.75,nangle=48)[
        "forward_pair_reaction_density_GeV4"]
    high=thermal_forward_rate_MB(T+h,61.75,nangle=48)[
        "forward_pair_reaction_density_GeV4"]
    energy_finite_difference=T*T*(high-low)/(2*h)
    analytic=equilibrium_pair_densities(T,nangle=48)[
        "event_total_energy_density_GeV5"]
    assert analytic>0
    assert abs(energy_finite_difference/analytic-1)<2e-4
