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
