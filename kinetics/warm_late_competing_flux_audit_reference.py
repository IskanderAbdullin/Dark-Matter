"""Audit frozen channel magnitudes and matched expansion endpoints."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from joint_visible_heavy_adiabatic_reference import snapshot_functions
from three_species_decay_onsager_reference import m,mp
from warm_corrected_late_prompt_cut_aligned_reference import excess
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom


def exact_free_W(s,b,C):
    """Liouville shift evaluated on the original state without refitting."""
    p=np.array(s["lab_momentum_edges_GeV"])
    B,_,F,_=snapshot_functions(s)
    x,w=roots_legendre(24)
    acc=0.
    for mass,deg,g in ((mp,1,B),(m,4,F)):
        cut=b*np.sqrt((C-mass)*(C+mass))
        edges=np.unique(np.r_[max(p[0],cut),p[p>cut]])
        dp=np.diff(edges)[:,None]
        u=edges[:-1,None]+dp*(x+1.)/2.
        W=np.maximum(0.,np.hypot(u/b,mass)-C)
        acc+=deg/(4.*np.pi**2*b**3)*np.sum(w*dp*u**2*g(np.hypot(u,mass))*W)
    return float(acc)


def audit(start, paths):
    r=read(start)
    s=r["snapshot"]
    assert r["hash"]==digest(s)
    e=[read(paths/path) for path in ("elastic-low.json","elastic-high.json")]
    p=[read(paths/path) for path in ("pair-low.json","pair-high.json")]
    assert all(len(t["rows"])==4 for t in e+p)
    assert len({(t["H_GeV"],t["normalization_rho_GeV4"],t["cut_GeV"]) for t in e+p})==1
    assert all(len({i["input_sha256"] for i in rows})==1 for rows in zip(*(t["rows"] for t in e+p)))
    direct=[]
    for rows in zip(*(t["rows"] for t in e+p)):
        f=rows[0]["prompt_signed_over_rhoH"]
        fid=rows[0]["input_sha256"]
        el=[i["elastic"][-1]["rates_over_rhoH"] for i in rows[:2]]
        pp=[i["pair_rates_over_rhoH"] for i in rows[2:]]
        largest_el=max(i[0]+i[1] for i in el)
        largest_pair=max(i[0]+i[1] for i in pp)
        assert all(i["elastic"][0]["rates_over_rhoH"]==[0.,0.,0.] for i in rows[:2])
        assert max(i[3] for i in pp)<3e-28
        assert (largest_el+largest_pair)/abs(f)<5e-6
        direct.append(dict(time=rows[0]["time"],input_sha256=fid,prompt=f,elastic_gross_and_signed=el,
                           pair_gross_signed_and_outgoing=pp,largest_sampled_gross_sum_relative_to_prompt=(
                           largest_el+largest_pair)/abs(f)))
    choices=("expanding-prompt.json","expanding-pair.json","expanding-all.json")
    dyn=[read(paths/x) for x in choices]
    assert [x["mode"] for x in dyn]==["prompt","prompt-pair","all"]
    assert all(x["input_hash"]==r["hash"] and digest(x["snapshot"])==x["hash"] for x in dyn)
    assert len({x["guide_hash"] for x in dyn})==1
    assert len({x["relative_log_expansion"] for x in dyn})==1
    C=dyn[0]["fixed_energy_cut_GeV"]
    assert C==e[0]["cut_GeV"] and all(x["fixed_energy_cut_GeV"]==C for x in dyn)
    h=dyn[0]["relative_log_expansion"]
    p0=np.asarray(s["lab_momentum_edges_GeV"])
    rho=get_mom(p0,[s["entropy_variable_"+x] for x in ("scalar","fermion")])[2].sum()
    init=exact_free_W(s,1.,C)/rho
    first=exact_free_W(s,np.exp(h/2),C)/rho
    first_error=max(abs(x["stages"][0]["W_over_initial_rho"]/first-1.) for x in dyn)
    assert first_error < 1e-10
    W=[]
    for x in dyn:
        f=x["snapshot"]
        assert np.array_equal(p0/np.exp(h),f["lab_momentum_edges_GeV"])
        z=sum(excess(f,C))/rho
        assert abs(x["stages"][-1]["W_over_initial_rho"]/z-1.)<2e-14
        assert x["solver"]["residual"]<8e-12
        W.append(z)
    return dict(frozen=direct,relative_log_expansion=h,exact_start_W_over_initial_rho=init,
                direct_first_redshift_W_over_initial_rho=first,first_reconstruction_relative_error=first_error,
                exact_free_endpoint_W_over_initial_rho=exact_free_W(s,np.exp(h),C)/rho,
                matched_endpoints=[dict(W_over_initial_rho=w,
                      signed_change_per_delta_log_a_over_initial_rho=(w-init)/h,
                      relative_difference_to_expanding_prompt=w/W[0]-1.,
                      mode=x["mode"],solver=x["solver"],state_hash=x["hash"],
                      pair_count_increment_over_start=x["pair_count_increment_over_start"],
                      relative_total_count_residual=x["relative_total_count_residual"],
                      relative_energy_residual=x["relative_energy_residual"])
                      for w,x in zip(W,dyn,strict=True)])


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("start","paths"):
        p.add_argument(key,type=Path)
    a=p.parse_args()
    print(json.dumps(audit(a.start,a.paths),indent=2,allow_nan=False))
