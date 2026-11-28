"""Combine disjoint warm pion sources without claiming a complete spectrum."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.special import gammaincc

from kloe_pion_band_portal_reference import POINT_MUON,table
from timelike_portal_onsager_coverage_reference import mpi
from warm_extra_injection_headroom_reference import run_case, pion_positive_Q
from warm_pion_bose_and_threshold_reference import low_strip_per_unit,bose


def read(p):
    return json.loads(Path(p).read_bytes())


def run(eos,kloe,kernel,babar,neutral,theory):
    K=read(kernel)["simultaneous_energy_density_coefficient_GeVminus2"]
    B=next(x for x in read(babar)["moment_variants"]
           if x["data"]=="BABAR 2026" and x["linear_error_multiple"]==2)
    N=read(neutral)["windows"][-1]["controls"][-1]["energy_per_unit_form_factor_squared"][1]
    Eb=B["gap_energy_upper_GeV4"]
    theory_amount=read(theory)["variants"][-1]["threshold_slice_energy_upper_GeV4"]
    s,cs,e=table(kloe)
    rmax=max((cs+e)*(s+.005)/POINT_MUON)
    assert rmax>0.
    u=.30**2
    Rstrip=.25*(1-4*mpi**2/u)**1.5
    B0=bose(.05)**2
    heat_unit=B0*gammaincc(5,mpi/.05)/6
    rows=[]
    for label,below,space,hat,theory_E in (
        ("measured charged pion and BABAR gap only",0.,0.,0.,0.),
        ("unit conditional pion form-factor caps",1.,1.,1.,0.),
        ("wider charged and neutral form-factor caps",4.,4.,10.,0.),
        ("dispersive threshold and wider elastic and neutral caps",0.,4.,10.,theory_amount)):
        controls=[]
        for n,du,ang in ((660,.012,8),(1200,.006,12)):
            def rates(T,data,splits):
                return (pion_positive_Q(T,data,splits)+below*low_strip_per_unit(T,u))*bose(T)**2
            c=run_case(eos,K,kloe,n,du,ang,pion_rate_override=rates,
                       pion_tail_multiplier=B0*(1.+below*Rstrip/rmax),
                       heat_multiplier=1.+space*heat_unit,return_components=True)
            C=c["baseline_components"]
            base=C["direct_Omega"]+C["recycle_Omega"]
            ecap=c["extra_energy_after_including_pions_to_touch_half_GeV4"]
            A,L=C["additional_early_portal_quadratic_coefficient"],C["additional_early_portal_linear_coefficient"]
            def upper(deposit):
                return base+L*deposit+A*deposit**2
            assert abs(upper(ecap)-.060)<3e-14
            eN=hat*hat*N
            assert below*theory_E == 0.  # Two alternatives for the same threshold strip.
            combined=Eb+eN+theory_E
            assert 0<combined<ecap
            controls.append(dict(time_panels=c["panels"],kernel_step=du,
                                 charged_including_below_without_BABAR_base=base,
                                 neutral_energy_GeV4=eN,
                                 BABAR_early_energy_upper_GeV4=Eb,
                                 dispersive_theory_early_energy_upper_GeV4=theory_E,
                                 joint_upper_scalar_Omega_h2=upper(combined),
                                 remaining_early_input_allowance_GeV4=ecap-combined,
                                 total_neutral_amplitude_budget_at_fixed_charged_caps=float(np.sqrt((ecap-Eb-theory_E)/N))))
        assert controls[-1]["joint_upper_scalar_Omega_h2"]<controls[0]["joint_upper_scalar_Omega_h2"]
        rows.append(dict(case=label,sub_point30_timelike_Fpi_squared_cap=below,
                         spacelike_Fpi_squared_cap=space,
                         neutral_Fhat_cap_GeVminus1=hat,controls=controls))
    return dict(scope="Initially empty 50 MeV, tabulated EOS. Distinct charged pion regions "
                      "below .30 (conditioned cap), .30–sqrt(.10) (whole BABAR moment cap), "
                      "and sqrt(.10)–sqrt(.95) (KLOE effective band). Extra neutral window and "
                      "elastic heating use separately stated form-factor caps. No full thermal "
                      "hadronic spectrum or abundance prediction.",
                input_sha256=[hashlib.sha256(Path(p).read_bytes()).hexdigest()
                              for p in (eos,kloe,kernel,babar,neutral,theory)],cases=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","kloe","kernel","babar","neutral","theory"):
        p.add_argument(k)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.kloe,a.kernel,a.babar,a.neutral,a.theory),indent=2,allow_nan=False))
