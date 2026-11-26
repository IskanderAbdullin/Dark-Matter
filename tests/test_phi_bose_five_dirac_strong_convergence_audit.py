"""First finite-time and spectral grid convergence for real five-field Yukawa."""
import math
import numpy as np
from kinetics.phi_bose_five_dirac_strong_convergence_audit import (
    construct,integrate_horizon,
)


def test_five_field_finitedt_first_order_on_multiple_physical_phi_energy_grids():
    for nphi in (8,12,16):
        m=construct(nphi,2.)
        initial=m.source(m.initial_etas())
        assert initial["H"]>0
        assert abs(initial["exact_event_U1"])<1e-18
        assert abs(initial["exact_event_energy"])<1e-17
        runs=[integrate_horizon(m,20000.,n) for n in (2,4,8)]
        n=[r["delta_Nphi_GeV3"] for r in runs]
        ratio=(n[0]-n[1])/(n[1]-n[2])
        print("PHYSICAL_5FIELD_TIME_CONVERGENCE",repr(dict(
            cutoff_GeV=2.,Nphi=nphi,dt_nominal_GeVinv=[10000.,5000.,2500.],
            delta_Nphi_GeV3=n,observed_first_order_ratio=ratio,
            entropy_changes=[x["delta_entropy_GeV3"] for x in runs],
            energy_defects=[x["absolute_delta_total_energy"] for x in runs],
            U1_defects=[x["absolute_U1_defect"] for x in runs])),flush=True)
        assert 1.4<ratio<2.6,(nphi,n,ratio)
        for v in runs:
            assert not v["rejected"]
            assert abs(v["absolute_U1_defect"])<3e-11
            assert abs(v["absolute_delta_total_energy"])<3e-11
            assert v["delta_entropy_GeV3"]>0
            assert v["min_step_entropy_certificate"]>-3e-10
            assert v["min_Bose_eta"]>0


def test_five_field_finite_time_not_a_spurious_source_of_net_U1_or_energy_at_P8():
    m=construct(8,8.)
    x=integrate_horizon(m,20000.,4)
    assert abs(x["absolute_U1_defect"])<3e-11
    assert abs(x["absolute_delta_total_energy"])<3e-11
    assert x["min_step_entropy_gain"]>0
    assert x["min_Bose_eta"]>0


def test_full_fivefield_physical_parent_UV_cutoff_sensitivity_at_same_finite_time():
    # Three hard parent cutoffs are independent PHYSICAL model truncations.
    # Strong-step positivity or tiny timestep residual is not evidence of
    # convergence of the full thermal parent UV tail.
    rows=[]
    for P in (2.,8.,16.):
        m=construct(12,P)
        src=m.source(m.initial_etas())
        out=integrate_horizon(m,20000.,4)
        rows.append(dict(P=P,deltaNphi=out["delta_Nphi_GeV3"],
                         rate_H=src["H"],energy_defect=out["absolute_delta_total_energy"],
                         U1_defect=out["absolute_U1_defect"],
                         entropy_change=out["delta_entropy_GeV3"],events=m.total_events))
        assert src["H"]>=0
        assert out["min_Bose_eta"]>0
        assert abs(out["absolute_U1_defect"])<3e-11
        assert abs(out["absolute_delta_total_energy"])<3e-11
    print("PHYSICAL_5FIELD_UV_CUTOFF",repr(rows),flush=True)
    assert all(v["entropy_change"]>0 for v in rows)


def test_full_pair_radial_Gauss_source_underresolution_at_large_hard_parent_cutoff():
    # This is a *physical source convergence* test, not a claim that 4
    # Gauss parent radial nodes at P=16 GeV sample thermal T~1 GeV.
    # Change only pair/decay radial counts; keep relative & azimuth/CM
    # fixed to separate radial from 2->2 scattering angular errors.
    rows=[]
    for nr in (4,6,8):
        m=construct(12,16.,nrad=nr,nrelative=8,ncm=32,naz=4)
        src=m.source(m.initial_etas())
        row=dict(P_GeV=16.,nrad=nr,events=m.total_events,
                 H=src["H"],rates={c["name"]:c["event_rate"] for c in src["channels"]},
                 U1=src["exact_event_U1"],E=src["exact_event_energy"])
        rows.append(row)
        assert src["H"]>0
        assert abs(src["exact_event_U1"])<1e-18
        assert abs(src["exact_event_energy"])<1e-17
    print("PHYSICAL_5FIELD_RADIAL_QUADRATURE",repr(rows),flush=True)
    assert all(v["events"]>0 for v in rows)


def test_independent_full_quantum_pair_radial_source_refines_beyond_n8():
    from kinetics.phi_bose_five_dirac_full_yukawa_weak import FullDiracYukawaWeak
    rows=[]
    for nr in (8,12,16):
        model=FullDiracYukawaWeak(
            P=16.,Ppair=16.,nphi=12,nrad=8,nangle=8,
            pair_nrad=nr,pair_nrelative=8,pair_ncm=32,pair_naz=4)
        res=model.source(T=(.9,1.05,1.,1.1,.98))
        ch=res["event_channels"]
        rows.append(dict(nrad=nr,
                         chi1_pair_rate=ch["chi1_pair"]["event_rate_GeV4"],
                         chi2_pair_rate=ch["chi2_pair"]["event_rate_GeV4"],
                         H=res["H_GeV4"],
                         events=model.event_counts))
        assert abs(res["vector_U1_defect_GeV4"])<1e-18
        assert abs(res["total_energy_defect_GeV5"])<1e-17
        assert res["H_GeV4"]>0
    print("PHYSICAL_PAIR_RADIAL_8_12_16",repr(rows),flush=True)
    # Positivity and conservation never imply convergence of a tiny
    # pair contribution hidden under the much stronger 1->2 decay H.
    assert all(v["chi1_pair_rate"]<0 for v in rows)
    assert all(v["chi2_pair_rate"]>0 for v in rows)


def test_fixed_positive_graded_pair_radial_source_convergence_and_exact_capacity():
    from kinetics.phi_bose_five_dirac_full_yukawa_weak import FullDiracYukawaWeak
    # Explicit physical geometry breakpoints, independent of occupations.
    P=16.
    edges=[0.,.25,1.,4.,P]
    rows=[]
    for nper in (2,3,4):
        m=FullDiracYukawaWeak(
            P=P,Ppair=P,nphi=12,nrad=8,nangle=8,
            pair_nrad=nper,pair_nrelative=8,pair_ncm=32,pair_naz=4,
            pair_radial_breaks=edges)
        assert m.pair1.radial_segments==m.pair2.radial_segments==4
        for pair in (m.pair1,m.pair2):
            assert len(pair.radial_p)==4*nper
            assert np.min(pair.radial_weights)>0
            capacity=float(np.dot(pair.radial_weights,pair.radial_p**2))
            assert math.isclose(capacity,P**3/3,rel_tol=3e-14)
        result=m.source(T=(.9,1.05,1.,1.1,.98))
        chans=result["event_channels"]
        assert abs(result["vector_U1_defect_GeV4"])<1e-18
        assert abs(result["total_energy_defect_GeV5"])<1e-17
        assert result["H_GeV4"]>0
        rows.append(dict(graded_nper=nper,actual_radial_nodes=4*nper,
                         chi1_pair_rate=chans["chi1_pair"]["event_rate_GeV4"],
                         chi2_pair_rate=chans["chi2_pair"]["event_rate_GeV4"],
                         total_H=result["H_GeV4"]))
    print("PHYSICAL_PAIR_GRADED_RADIAL_NODES",repr(rows),flush=True)


def test_graded_pair_radial_integrator_keeps_full_fivefield_finite_time_U1_E_entropy():
    from kinetics.phi_bose_five_dirac_full_yukawa_implicit import FullDiracPhysicalMoments
    from kinetics.phi_bose_five_dirac_kl_certificate import quantum_KL_certificate
    m=FullDiracPhysicalMoments(
        P=16.,Ppair=16.,N1=3,N2=3,Nphi=12,nq=18,
        nrad=4,nangle=8,pair_nrad=2,
        pair_nrelative=8,pair_ncm=32,pair_naz=4,
        pair_radial_breaks=[0.,.25,1.,4.,16.])
    old=m.initial_etas()
    result=m.adaptive_step(old,5000.)
    assert result["min_Bose_eta"]>0
    assert result["entropy_increase"]>0
    assert result["entropy_certificate"]>=-1e-10
    assert abs(result["energy_defect"])<3e-11
    assert abs(result["U1_defect"])<3e-11
    cert=quantum_KL_certificate(m,old,result["etas"],result["accepted_dt"])
    assert abs(cert["exact_entropy_decomposition_defect"])<3e-10
