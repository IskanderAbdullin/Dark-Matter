"""Five-dimensional arbitrary-isotropic nonthermal QED weak-moment tests."""
import math
import numpy as np
from scipy.special import expit
from kinetics.qed_pair_nonthermal_weak import weak_moments
from kinetics.qed_pair_full5d import thermal_collision

def test_independent_thermal_reference_nonzero_affinity():
    T=.1
    f=lambda E:expit((.01*T-E)/T)
    b=lambda E:1/np.expm1(E/T)
    a=weak_moments(f,f,b,Escale=T,nrad=24,nang=16)
    ref=thermal_collision(T,T,mu_minus=.01*T,mu_plus=.01*T,
                          nrad=24,nang=16)
    assert math.isclose(a["R_net_GeV4"],ref["R_net_GeV4"],
                        rel_tol=2e-7)
    assert math.isclose(a["entropy_production_GeV4"],
                        ref["entropy_production_GeV4"],rel_tol=2e-7)

def test_thermal_fixed_point_has_negligible_net_and_entropy():
    T=.1;f=lambda E:expit(-E/T);b=lambda E:1/np.expm1(E/T)
    a=weak_moments(f,f,b,Escale=T,nrad=24,nang=16)
    assert abs(a["R_net_GeV4"])<1e-14*a["R_forward_GeV4"]
    assert 0<=a["entropy_production_GeV4"]<1e-22*a["R_forward_GeV4"]

def test_nonthermal_electron_packet_gain_loss_and_h_theorem():
    T=.1
    f0=lambda E:expit(-E/T)
    for strength,sign in [(.15,1),(-.15,-1)]:
        fminus=lambda E:expit(-E/T+strength*np.exp(-E/T))
        b=lambda E:1/np.expm1(E/T)
        a=weak_moments(fminus,f0,b,Escale=T,nrad=24,nang=16)
        assert sign*a["R_net_GeV4"]>0
        assert a["entropy_production_GeV4"]>0
        assert abs(a["Q_total_residual_GeV5"])<1e-12*abs(a["Q_photon_GeV5"])
        assert a["N_photon_GeV4"]==-2*a["N_e_minus_GeV4"]

def test_nonthermal_photon_peak_can_reverse_net_while_entropy_positive():
    T=.1
    f=lambda E:expit(-E/T)
    thermal_b=lambda E:1/np.expm1(E/T)
    for amplitude,sign in [(+.2,-1),(-.2,+1)]:
        b=lambda E:thermal_b(E)*(1+amplitude*np.exp(-((E-1.1)/.2)**2))
        a=weak_moments(f,f,b,Escale=T,nrad=24,nang=16)
        assert sign*a["R_net_GeV4"]>0
        assert a["entropy_production_GeV4"]>0
        assert abs(a["Q_total_residual_GeV5"])<1e-12*abs(a["Q_photon_GeV5"])
