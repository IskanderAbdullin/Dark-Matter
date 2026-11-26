"""Independent 5D Pauli/Bose QED thermal reaction and conservation tests."""
import math
from kinetics.qed_pair_full5d import thermal_collision
from kinetics.qed_pair_annihilation_control import initial_fd_loss_density

def test_vacuum_5d_amplitude_loss_reproduces_independent_total_sigma():
    for T in (.025,.1):
        a=thermal_collision(T,None,nrad=32,nang=24,nphi=12)
        ref=initial_fd_loss_density(T,n=48,cut=32)
        assert math.isclose(a["R_forward_GeV4"],ref,rel_tol=2e-8)
        assert a["R_inverse_GeV4"]==0
        assert abs(a["energy_residual_GeV5"])<1e-12*abs(a["Q_photon_GeV5"])
        assert a["electric_charge_source_GeV4"]==0

def test_thermal_pair_chemical_equilibrium_even_with_nonzero_charge_mu():
    T=.1
    a=thermal_collision(T,T,mu_minus=.02,mu_plus=-.02,nrad=20,nang=16)
    assert a["R_net_GeV4"]==0.
    assert a["entropy_production_GeV4"]==0.
    assert math.isclose(a["R_forward_GeV4"],a["R_inverse_GeV4"],
                        rel_tol=5e-13)

def test_chemical_affinity_entropy_exact_and_conservation():
    T=.1; A_over_T=.02
    a=thermal_collision(T,T,mu_minus=.001,mu_plus=.001,
                        nrad=24,nang=20)
    assert a["R_net_GeV4"]>0
    assert math.isclose(a["R_net_GeV4"]/a["R_forward_GeV4"],
                        -math.expm1(-A_over_T),rel_tol=1e-12)
    assert math.isclose(a["entropy_production_GeV4"],
                        A_over_T*a["R_net_GeV4"],rel_tol=1e-12)
    assert abs(a["energy_residual_GeV5"])<1e-12*abs(a["Q_photon_GeV5"])
    assert a["N_e_minus_source_GeV4"]==a["N_e_plus_source_GeV4"]
    assert a["N_photon_source_GeV4"]==-2*a["N_e_minus_source_GeV4"]

def test_hotter_photons_create_pairs_cold_photons_annihilate_pairs():
    cold=thermal_collision(.1,.09,nrad=20,nang=16)
    hot=thermal_collision(.1,.11,nrad=20,nang=16)
    assert cold["R_net_GeV4"]>0 and hot["R_net_GeV4"]<0
    assert cold["entropy_production_GeV4"]>0 and hot["entropy_production_GeV4"]>0
    for a in (cold,hot):
        assert abs(a["energy_residual_GeV5"])<2e-12*abs(a["Q_photon_GeV5"])

def test_radial_refinement_of_full_5d_chemical_source():
    T=.1
    a=thermal_collision(T,T,mu_minus=.001,mu_plus=.001,
                        nrad=24,nang=20)
    b=thermal_collision(T,T,mu_minus=.001,mu_plus=.001,
                        nrad=32,nang=20)
    assert abs(a["R_net_GeV4"]/b["R_net_GeV4"]-1)<5e-8
