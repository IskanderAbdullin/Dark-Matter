"""Conditional per-spin overshoot sensitivity at the archived late q64 point."""

import argparse
import json
from pathlib import Path

from scipy.optimize import brentq

from three_species_decay_onsager_reference import M,m
from warm_isolated_recycling_energy_envelope_reference import radiation_bound
from warm_light_occupation_stage_audit_reference import read
from warm_signed_fermion_piecewise_clock_reference import weights
from warm_signed_fermion_time_comparison_reference import setup
from warm_signed_quadratic_overshoot_reference import inverse,lower


def run(eos,D):
    a=read(D/"warm-signed-fermion-piecewise-clock-lower-2026-11-26.json")["cases"]
    Kpi=read(D/"warm-signed-fermion-with-charged-pions-2026-11-26.json")["cases"][-1]["controls"][-1]["pion_pair_sigma_v_cap_with_stated_support"]
    K=read(D/"warm-fermion-survival-with-leptons-2026-11-26.json")["coefficient"]
    J=sum(inverse(x,y,620,1240)/x**2 for x,y in ((m,M),(M,m)))
    answer=[]
    for S in a:
        T=S["T_GeV"]; N=S["number_at_anchor_GeV3"]
        E=S["energy_upper_GeV4"]; z=S["unconverted_start_proxy"]
        f=radiation_bound(eos,T,500.)
        _,B,*r=setup(N,E,K,f,2800)
        W=weights(eos,T,500.,2800,f["photon_anchor_GeV"])
        def objective(d):
            return z*lower((W,B+.5*N*Kpi,*r),d,E*E*J/N)
        assert abs(objective(.5)/S["controls"][-1]["lower_Omega_h2"][-1]["local_floor_lower"]-1)<1e-13
        root=brentq(lambda d:objective(d)-.12,.5,.51,xtol=1e-12)
        samples=[dict(cap=d,lower=objective(d)) for d in (.5,.5001,.5003,.5005,.5006)]
        answer.append(dict(T_GeV=T,last_cap_certifying_exclusion_approximately=root,sample=samples))
    return dict(scope="Same selected vacuum network, scalar not capped. Bounds the extra potentially "
                      "positive Bose-quadratic loss for small postulated fermion cap overshoot. "
                      "No dynamical cap or robustness up to full Pauli saturation proved.",cases=answer)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("eos",type=Path)
    p.add_argument("--tapes",type=Path,default=Path("results"));a=p.parse_args()
    print(json.dumps(run(a.eos,a.tapes),indent=2,allow_nan=False))
