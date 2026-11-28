"""An exclusion-envelope map, not a fit, for two fixed warm source models."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from warm_pion_bose_and_threshold_reference import run as pion_fold


def upper(D,R,J,epsilon_over_ref,y_over_ref):
    x=epsilon_over_ref**2
    return D*x + R*x*x*y_over_ref**4*np.exp(2*J*(x-1))


def study(name,D,R,J,need=.060):
    # For the stated fixed-width model, these curves are monotone in each
    # coupling norm; no physical root or rate equality is asserted.
    answer=[]
    for e in (.7,.85,1.,1.1,1.2,1.35,1.5):
        assert upper(D,R,J,e,0.) < need
        f=brentq(lambda y: upper(D,R,J,e,y)-need,0.,5.,xtol=1e-13)
        assert abs(upper(D,R,J,e,f)-need)<1e-12
        x=e*e
        analytic=((need-D*x)/(R*x*x*np.exp(2*J*(x-1))))**.25
        assert abs(f/analytic-1.)<2e-13
        answer.append(dict(epsilon_over_reference=e,y_over_reference_for_upper_to_touch_half=f,
                           upper_at_unmodified_y=float(upper(D,R,J,e,1.))))
    e=brentq(lambda e: upper(D,R,J,e,1.)-need,0.,5.,xtol=1e-13)
    return dict(model=name,original_direct=D,original_recycling=R,original_heat_exponent=J,
                original_upper=D+R,epsilon_over_reference_for_upper_to_touch_half_at_fixed_y=e,
                y_over_reference_for_upper_to_touch_half_at_fixed_epsilon=answer[2]["y_over_reference_for_upper_to_touch_half"],
                grid=answer)


def run(eos,pion_data,kernel,baseline,prior_pion):
    z=json.loads(Path(baseline).read_bytes())["tabulated_EOS_controls"][-1]
    p=pion_fold(eos,pion_data,kernel,only="Bose_and_widened_threshold_and_elastic_caps",
                return_components=True)["cases"][0]["comparison"]
    pc=p["baseline_components"]
    old=json.loads(Path(prior_pion).read_bytes())["cases"][-1]["comparison"]
    assert abs(pc["direct_Omega"]+pc["recycle_Omega"]-old[
        "phased_measured_pion_and_lepton_scalar_Omega_ceiling"])<2e-15
    return dict(scope="Conditional instantaneous 50 MeV vacuum benchmark with fixed "
                      "mediator width and visible clock. Changes only the source "
                      "coupling and the specified dark Yukawa kernel, not other reactions.",
                epsi_ref=1e-5,y_ref=1e-3,reference_half=.060,
                input_sha256=[hashlib.sha256(Path(p).read_bytes()).hexdigest()
                              for p in (eos,pion_data,kernel,baseline,prior_pion)],
                cases=[study("vacuum_leptons_only",z["Omega_phi_direct_count_ceiling"],
                             z["Omega_phi_recycle_ceiling"], z["heating_exponent"]),
                       study("with_selected_pions_and_stated_caps",pc["direct_Omega"],
                             pc["recycle_Omega"],pc["heat_exponent"])])


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","pion_data","kernel","baseline","prior_pion"):
        p.add_argument(k)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.pion_data,a.kernel,a.baseline,a.prior_pion),indent=2,allow_nan=False))
