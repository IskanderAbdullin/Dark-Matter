"""Optional positive geometric physical phi first-cell DOS quadrature
integrated into the REAL five-population on-shell 4-channel quantum solver.

These tests target actual physical entropy moments and G, not simply weak
collision rates, and do not establish BEC/SM/FRW/relic.
"""
import math
import numpy as np
from kinetics.phi_bose_five_dirac_full_yukawa_implicit import FullDiracPhysicalMoments
from kinetics.phi_bose_five_dirac_kl_certificate import quantum_KL_certificate
from kinetics.phi_bose_physical_cutcell_weak import MPHI,M1,M2
from kinetics.phi_bose_critical_ir_susceptibility import stretched_physical_moments


def _model(IR):
    return FullDiracPhysicalMoments(
        P=2.,Ppair=2.,N1=3,N2=3,Nphi=8,nq=18,
        nrad=4,nangle=8,pair_nrad=4,pair_nrelative=8,
        pair_ncm=32,pair_naz=4,phi_IR_composite=IR)


def _exact_gap_etas(model,delta=1e-16):
    eta=model.initial_etas(T=(1.,)*5)
    # Given directly in terms of positive physical gap, avoid roundoff
    # subtracting very close mphi and mu at the first Bose node.
    phi=model.species[4]
    eta[4]=(phi.nodal_energies-MPHI+delta)
    assert eta[4][0]==delta and min(eta[4])>0
    return eta


def test_first_physical_massive_Bose_cell_G_and_N_improve_without_changing_capacity():
    std=_model(False)
    fixed=_model(True)
    assert fixed.phi_IR_composite and fixed.phi_IR_composite_nodes==168
    assert not std.phi_IR_composite
    exact=stretched_physical_moments(1e-16)
    orig_phi=std.species[4]
    new_phi=fixed.species[4]
    eta_std=_exact_gap_etas(std)
    eta_new=_exact_gap_etas(fixed)
    raw=orig_phi.evaluate(eta_std[4],with_G=True)
    better=new_phi.evaluate(eta_new[4],with_G=True)
    G_old=float(raw["G"][0,0,0])
    G_new=float(better["G"][0,0,0])
    assert G_old/exact["susceptibility_G00_GeV3"]<.004
    assert abs(G_new/exact["susceptibility_G00_GeV3"]-1)<2e-7
    assert abs(better["M"][0]/exact["first_number_moment_M0_GeV3"]-1)<3e-10
    assert math.isclose(new_phi.edges[1],orig_phi.edges[1],rel_tol=2e-13)
    assert np.min(new_phi.weights[0])>0
    assert np.min(new_phi.H[0])>=0
    for sp in fixed.species:
        p0=np.sqrt(np.maximum(0,(sp.edges[:-1]-sp.mass)*(sp.edges[:-1]+sp.mass)))
        p1=np.sqrt(np.maximum(0,(sp.edges[1:]-sp.mass)*(sp.edges[1:]+sp.mass)))
        cap=sp.g*(p1**3-p0**3)/(6*math.pi**2)
        got=np.sum(sp.weights,axis=1)
        assert np.max(abs(got/cap-1))<5e-13


def test_optional_static_grid_does_not_change_any_actual_on_shell_1to2_or_2to2_event_rate():
    std=_model(False)
    fixed=_model(True)
    e0=_exact_gap_etas(std,1e-12)
    e1=_exact_gap_etas(fixed,1e-12)
    src0=std.source(e0)
    src1=fixed.source(e1)
    assert std.total_events==fixed.total_events
    for a,b in zip(src0["channels"],src1["channels"]):
        assert a["name"]==b["name"]
        assert math.isclose(a["event_rate"],b["event_rate"],
                            rel_tol=2e-12,abs_tol=1e-22)
    assert math.isclose(src0["H"],src1["H"],rel_tol=2e-12,abs_tol=1e-22)
    # First cell truly collision-inaccessible at hard P2 for all FOUR channels.
    for src in (src0,src1):
        assert max(abs(src["C"][fixed.offsets[4]:fixed.offsets[4]+3]))==0


def test_positive_composite_phi_quadrature_keeps_exact_KL_residual_work_identity():
    m=_model(True)
    old=_exact_gap_etas(m,1e-12)
    new=[e.copy() for e in old]
    new[0]+=0.003*np.sin(np.arange(len(new[0]))+.4)
    new[4][0]*=1.1
    cert=quantum_KL_certificate(m,old,new,1000.)
    assert all(v>=-5e-12 for v in cert["physical_quantum_KL_by_species"])
    assert abs(cert["exact_entropy_decomposition_defect"])<3e-10
    assert math.isfinite(cert["entropy_gain"])


def test_subcritical_nearcritical_phi_optional_physical_root_and_KL():
    m=_model(True)
    old=_exact_gap_etas(m,1e-16)
    trial=m.adaptive_step(old,5000.)
    assert 0<trial["accepted_dt"]<=5000.
    assert trial["min_Bose_eta"]>0
    assert trial["scaled_residual"]<2e-7
    # Critical full-S and full-E differ only below their float64 ULP.
    # Accepted means no *RESOLVED* entropy or conservation violation;
    # it does NOT imply a measurable positive entropy increment.
    assert trial["roundoff_entropy_floor"]>0
    assert trial["roundoff_energy_floor"]>0
    assert trial["roundoff_U1_floor"]>0
    assert not trial["entropy_certificate_sign_resolved"]
    assert not trial["entropy_gain_sign_resolved"]
    assert abs(trial["entropy_certificate"])<trial["roundoff_entropy_floor"]
    assert abs(trial["energy_defect"])<trial["roundoff_energy_floor"]
    assert trial["entropy_increase"]>=-1e-10
    assert trial["entropy_certificate"]>=-1e-10
    assert abs(trial["energy_defect"])<3e-11
    assert abs(trial["U1_defect"])<3e-11
    certificate=quantum_KL_certificate(m,old,trial["etas"],trial["accepted_dt"])
    assert abs(certificate["exact_entropy_decomposition_defect"])<3e-10
    assert certificate["microscopic_H_new"]>=0
    # Check the nonnegative quantum entropy-dissipation event measure
    # separately, even when direct total ΔS subtraction is inconclusive.
    assert all(x>=-1e-11 for x in
               certificate["physical_quantum_KL_by_species"])
