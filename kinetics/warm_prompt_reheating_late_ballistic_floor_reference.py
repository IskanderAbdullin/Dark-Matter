"""Conditional collisionless-redshifting loss-only floor from the two-temperature 2 MeV handoff."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import expit,roots_legendre

from archived_collisionless_streaming_reference import Background
from conditional_relic_normalization_reference import conversion
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import m
from warm_pair_shell_survival_bound_reference import exposure
from warm_prompt_reheating_expanding_slow_reference import visible


def run(eos,prior):
    d=json.loads(prior.read_text())
    x=next(z for z in d["controls"][-1]["heat_arms"] if z["prompt_heat"])["rows"][-1]
    b=Background(eos)
    tv=x["visible_T"]
    T=x["light_and_scalar_dark_T"][0]
    peak=x["light_and_scalar_phase_peak"][0]
    _,s=visible(b,tv)
    initial=x["yields"][0]*s
    Hph=np.sqrt(np.pi**2/45)*tv*tv/MPL
    cv=conversion(eos)["Omega_h2_per_GeV_times_yield"]*m
    phases=[]
    for n in (320,960):
        p=np.geomspace(T*1e-5,.035,n+1)
        p=np.r_[0.,p]
        K=lambda q:q*q/(np.hypot(q,m)+m)
        f=expit(np.log(peak)-K(p)/T)
        a,h=p[:-1],p[1:]
        weights=(h-a)*(h*h+h*a+a*a)*4/(6*np.pi**2)
        Nlo,Nhi=weights*f[1:],weights*f[:-1]
        gx,gw=roots_legendre(14)
        pq=(h+a)[:,None]/2+(h-a)[:,None]*gx/2
        N0=float(np.sum((h-a)[:,None]*gw*expit(np.log(peak)-K(pq)/T)*pq*pq)/np.pi**2)
        assert Nlo.sum()<N0<Nhi.sum() and abs(N0/initial-1)<1e-9
        entries=[]
        for cap in (.5,1.,4.,10.):
            exp=exposure(h,Nhi,Hph,cap)
            # Using shell upper populations in the loss subtraction
            # avoids creating an artificial normalization floor.
            lost=float(Nhi@(-np.expm1(-exp)))
            kept=max(0.,N0-lost)
            entries.append(dict(assumed_scalar_occupation_supremum=cap,retained_lower_fraction=kept/N0,
                                initial_omega_h2_if_passively_retained=N0/s*cv,
                                late_conditional_omega_h2_lower=kept/s*cv))
        phases.append(dict(momentum_panels=n,pmax_GeV=.035,controls=entries))
    return dict(scope="Charge-symmetric light fermions in freely redshifting momentum cells with loss only, "
                      "hard retained cutoff 35 MeV and no exterior or new population, elastic momentum transport "
                      "or signed inverse/source. The final scalar Bose factor is limited by an assumed uniform cap "
                      "on the whole later history. Neither this cap nor the cold phase history is shown to hold "
                      "in the actual model. The positive shellwise floor uses photons alone for the expansion bound.",
                visible_anchor_GeV=tv,fermion_temperature_GeV=T,controls=phases)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("prior",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.prior),indent=2,allow_nan=False))
