"""Reallocate a dispersive theory moment to the remaining sub-0.30 GeV range."""
import argparse
import hashlib
import json
from pathlib import Path

from timelike_portal_onsager_coverage_reference import mpi
from warm_babar_dispersive_gap_reference import thermal_window_upper


def run(eos,babar):
    prev=json.loads(Path(babar).read_bytes())
    E=next(z["gap_energy_upper_GeV4"] for z in prev["moment_variants"]
           if z["data"]=="BABAR 2026" and z["linear_error_multiple"]==2)
    ecap=prev["prior_remaining_after_measured_bosonic_band_GeV4"]
    base=prev["prior_scalar_Omega_upper"]
    smin,smax=(2*mpi)**2,.30**2
    folds=[thermal_window_upper(eos,n,t,smin,smax) for n,t in [(40,550),(90,1200)]]
    assert folds[-1]["Q_over_positive_moment_upper_GeV4"]<folds[0]["Q_over_positive_moment_upper_GeV4"]
    inputs=[]
    for extra in (0.,1.,2.):
        # The e+e- only theory determination, not a numerical data table.
        I=(132.91+.76*extra)*1e-10
        q=folds[-1]["Q_over_positive_moment_upper_GeV4"]*I
        share=(E+q)/ecap
        assert 0.<share<1.
        inputs.append(dict(prior_measurement_gap_E_upper_GeV4=E,
                           theory_error_multiple=extra,
                           theory_full_positive_moment=I,
                           threshold_slice_energy_upper_GeV4=q,
                           remaining_headroom_GeV4=ecap-E-q,
                           total_budget_share_spent=share,
                           scalar_Omega_secant_upper_with_both=base+share*(.060-base)))
    return dict(scope="Conditional vacuum theory prior on the positive pion moment below .63 GeV, "
                      "reallocated only to 2m_pi–.30 GeV. Separate BABAR prior is used only on "
                      ".30–sqrt(.10) GeV; KLOE is already above sqrt(.10). Not a new measurement "
                      "or a statistical combination, no full spectrum produced.",
                input_sha256=[hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in (eos,babar)],
                threshold_s_GeV2=[smin,smax],folds=folds,variants=inputs)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos")
    p.add_argument("babar")
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.babar),indent=2,allow_nan=False))
