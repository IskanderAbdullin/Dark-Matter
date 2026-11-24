"""Physical weak five-field full Dirac Yukawa: two conjugate 1↔2 + BOTH 2↔2.

These tests establish simultaneous real-mass N/E/U1/H weak bookkeeping,
not finite-time strong quantum, SM/FRW or actual frozen 86eV DM.
"""
import math
import numpy as np
from kinetics.phi_bose_five_dirac_full_yukawa_weak import (
    FullDiracYukawaWeak,CHARGES,
)
from kinetics.phi_bose_pair_physical_onshell_quantum import (
    differential_pair_sigma,PhysicalPairEvents,
)
from kinetics.phi_bose_physical_cutcell_weak import M1,M2,MPHI


def _m():
    return FullDiracYukawaWeak(P=2.,Ppair=2.,nphi=8,nrad=4,nangle=8,
                              pair_nrad=4,pair_nrelative=8,
                              pair_ncm=32,pair_naz=4)


def test_both_dirac_charge_conjugate_decays_reproduce_old_physical_3field_weak_normalization():
    m=_m()
    e=m.decay.initial_etas(T=(.9,1.,1.))
    old=m.decay.source(e)
    new=m.source(T=(.9,.9,1.,1.,1.),
                 include_decay=True,include_pairs=False)
    assert all(new["event_channels"][k]["event_count"]==len(m.decay.W)
               for k in ("chi1_decay","anti_chi1_decay"))
    assert math.isclose(new["dN"][4],old["Rphi"],rel_tol=3e-14)
    assert math.isclose(new["dQ"][4],old["Qphi"],rel_tol=3e-14)
    assert math.isclose(new["dQ"][0]+new["dQ"][1],old["Q1"],rel_tol=3e-14)
    assert math.isclose(new["dQ"][2]+new["dQ"][3],old["Q2"],rel_tol=3e-14)
    assert math.isclose(new["H_GeV4"],old["H"],rel_tol=3e-14)
    assert abs(new["vector_U1_defect_GeV4"])<1e-25
    assert abs(new["total_energy_defect_GeV5"])<5e-24


def test_full_five_species_common_T_U1_chemical_equilibrium_exact():
    m=_m()
    v=m.source(T=(1.,)*5,mu=(.05,-.05,.05,-.05,0.))
    assert set(v["event_channels"])=={
        "chi1_decay","anti_chi1_decay","chi1_pair","chi2_pair"}
    for c in v["event_channels"].values():
        assert c["min_weight"]>0
        assert abs(c["event_rate_GeV4"])<1e-24
        assert abs(c["quantum_H_GeV4"])<1e-26
    assert abs(v["H_GeV4"])<1e-26
    assert abs(v["vector_U1_defect_GeV4"])<1e-28
    assert abs(v["total_energy_defect_GeV5"])<1e-23


def test_swapped_heavy_Dirac_initial_pair_has_real_positive_onshell_2scalar_source():
    # Unlike the original light chi2 pair, this physically allowed channel
    # has initial mχ1=.1GeV and t/u virtual χ2=.02GeV.
    m=PhysicalPairEvents(P=2.,nrad=4,nrelative=8,ncm=32,naz=4,
                         m_in=M1,m_virtual=M2)
    assert min(m.W)>0
    assert m.max_kinematic_energy_defect<1e-13
    assert min(m.e_phi3)>MPHI and min(m.e_phi4)>MPHI
    assert m.MB_forward(1.)>0
    for E in (.21,.4,.8,1.5):
        u=np.array([-.8,0.,.5,.95])
        d=differential_pair_sigma(E*E,u,m_in=M1,m_virtual=M2)
        assert np.all(d>0)
    q=m.quantum_weak(muchi=.05,muanti=-.05,muphi=-.2)
    assert q["quantum_H_GeV4"]>0
    assert q["R_events_GeV4"]>0


def test_full_nonthermal_microevent_H_and_true_vector_charge_and_energy():
    m=_m()
    s=m.source(T=(.9,1.05,1.,1.1,.98),
               curvature=(.02,.01,.015,-.01,.005))
    assert s["H_GeV4"]>0
    assert all(ch["min_micro_JA"]>=-1e-11
               for ch in s["event_channels"].values())
    assert math.isclose(np.dot(CHARGES,s["dN"]),0.,abs_tol=5e-24)
    assert abs(s["vector_U1_defect_GeV4"])<5e-24
    assert abs(s["total_energy_defect_GeV5"])<5e-22
    assert s["class_ready"] is False


def test_old_particle_only_nonzero_mu_fixed_point_is_not_full_chemical_equilibrium():
    m=_m()
    s=m.source(T=(1.,)*5,mu=(-.15,.15,.05,-.05,-.2))
    assert abs(s["event_channels"]["chi1_decay"]["event_rate_GeV4"])<1e-25
    for name in ("anti_chi1_decay","chi1_pair","chi2_pair"):
        assert s["event_channels"][name]["event_rate_GeV4"]>0
        assert s["event_channels"][name]["quantum_H_GeV4"]>0
    assert s["H_GeV4"]>0


def test_instantaneous_full_network_is_not_fake_two_charge_closed_one_channel():
    m=_m()
    s=m.source(T=(1.,)*5,mu=(0.,0.,0.,0.,-.2),
               include_decay=False,include_pairs=True)
    # Pair number change breaks old N1+Nphi conservation.
    assert s["dN"][4]>0
    assert s["dN"][0]<0 and s["dN"][1]<0
    assert s["dN"][2]<0 and s["dN"][3]<0
    assert math.isclose(s["dN"][0],s["dN"][1],rel_tol=2e-14)
    assert math.isclose(s["dN"][2],s["dN"][3],rel_tol=2e-14)
    assert math.isclose(s["dN"][4],
                        -sum(s["dN"][:4]),rel_tol=2e-14)
    # The net vector U1, rather than the positive particle counts, survives.
    assert abs(s["vector_U1_defect_GeV4"])<1e-24
