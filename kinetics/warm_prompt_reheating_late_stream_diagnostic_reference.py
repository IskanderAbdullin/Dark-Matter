"""Compare the 2 MeV conditional thermal shape and its rest-density proxy."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from archived_collisionless_streaming_reference import Background
from conditional_relic_normalization_reference import conversion
from three_species_decay_onsager_reference import m,mp


def run(eos,prior):
    bg=Background(eos)
    d=json.loads(prior.read_text())
    x=next(p for p in d["controls"][-1]["heat_arms"] if p["prompt_heat"])["rows"][-1]
    Td=x["light_and_scalar_dark_T"][1]
    tv=x["visible_T"]
    peak=x["light_and_scalar_phase_peak"][1]
    cv=conversion(eos)["Omega_h2_per_GeV_times_yield"]
    out=[]
    for n in (12,22):
        cuts=np.r_[0.,np.geomspace(1e-6,72.,38)]
        a,b=(cuts[:-1]+cuts[1:])/2,np.diff(cuts)/2
        r,w=roots_legendre(n)
        p=Td*(a[:,None]+b[:,None]*r)
        K=p*p/(np.hypot(p,mp)+mp)
        F=peak*np.exp(-K/Td)
        mom=Td*w*b[:,None]*p*p*F/(1-F)
        q=p.ravel()*bg.scale(tv)
        weight=mom.ravel()/mom.sum()
        dir_=bg.streams(q,mp,bg.scale(tv),bg.aeq,8)
        out.append(dict(radial=n,rms_q_over_Tgamma0=float(np.sqrt(weight@(q*q))/bg.T0),
                        mean_q_over_Tgamma0=float(weight@q/bg.T0),
                        average_path_after_2MeV_until_equality_Mpc=float(weight@dir_),
                        fraction_p_over_Tgamma_at_2MeV_below_point_one=float(weight@(p.ravel()/tv<.1))))
    assert abs(out[0]["average_path_after_2MeV_until_equality_Mpc"]/out[-1]["average_path_after_2MeV_until_equality_Mpc"]-1)<2e-8
    om=np.array(x["yields"])*np.array([m,mp])*cv
    return dict(scope="Normal-phase scalar shape from the isolated two-T probe at 2 MeV. The ballistic path uses a fixed "
                      "reference background and assumes particles free-stream after the handoff, applicable as a shape "
                      "diagnostic only: the same run has a large surviving stable rest mass and residual annihilation. "
                      "No self-consistent expansion, final yield, transfer or data comparison.",
                Omega_h2_light_scalar_if_counts_hypothetically_freeze_at_anchor=om.tolist(),
                scalar_mass_fraction_of_both_frozen_proxies=om[1]/om.sum(),controls=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("prior",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.prior),indent=2,allow_nan=False))
