"""First finite-time and spectral grid convergence for real five-field Yukawa."""
import math
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
