"""One-way on-shell exterior birth rates on current archived warm shapes."""

import argparse
import json
from pathlib import Path

import numpy as np

from late_prompt_exterior_snapshot_reference import outward
from scalar_pair_expanding_continuation_reference import moments,e_total
from three_species_decay_onsager_reference import m,mp,width
from warm_light_occupation_stage_audit_reference import digest,read


def crop(s,edge):
    p=np.asarray(s["lab_momentum_edges_GeV"])
    j=int(np.argmin(abs(p-edge)))
    assert abs(p[j]/edge-1)<1e-12 and 0<j<len(p)-1
    return dict(s,lab_momentum_edges_GeV=p[:j+1].tolist(),
                entropy_variable_scalar=np.asarray(s["entropy_variable_scalar"])[:,:j].tolist(),
                entropy_variable_fermion=np.asarray(s["entropy_variable_fermion"])[:,:j].tolist())


def normed(S,H,settings):
    U=moments(np.asarray(S["lab_momentum_edges_GeV"]),
            [S["entropy_variable_scalar"],S["entropy_variable_fermion"]])
    N=U[:,0].sum(axis=1)
    rho=e_total(U)
    rows=[]
    for o in settings:
        a=outward(S,1.,o)
        assert np.all(a[:6]>0.) and np.all(a[6:]==0.)
        total_GeV4=(a[4:6]+a[2:4]*np.array([mp,m]))/H
        rows.append(dict(order=o,outer_positive_number_rate_GeV4=a[2:4].tolist(),
                         fraction_of_gross_decay_by_species=(a[2:4]/a[1]).tolist(),
                         rate_over_current_species_count_per_Hubble=(a[2:4]/(N*H)).tolist(),
                         outer_full_energy_over_retained_dark_energy_per_Hubble=(total_GeV4/rho).tolist(),
                         gross_decay_relative_difference_from_unblocked_inverse=(a[1]/a[0]-1.)))
    keys=("outer_positive_number_rate_GeV4","outer_full_energy_over_retained_dark_energy_per_Hubble")
    assert all(np.max(abs(np.asarray(rows[-1][k])/rows[-2][k]-1))<5e-5 for k in keys)
    return rows


def run(D):
    ref=read(D/"warm-near-rest-light-inverse-decay-2026-11-26.json.gz")
    edges=[]; rows=[]
    for q in (16,32):
        s=read(D/f"warm-light-occupation-next-q{q}-tape-2026-11-26.json.gz")
        s+=[z["checkpoint"] for z in read(D/f"warm-observed-suboperator-q{q}-tape-2026-11-26.json.gz")["blocks"]]
        match=next(z for z in ref["cases"] if z["q"]==q)["stages"]
        result=[]
        for j,(r,d) in enumerate(zip(s,match,strict=True)):
            S=r["snapshot"]
            assert r["snapshot_sha256"]==d["snapshot_sha256"]==digest(S)
            H=d["fine"]["H_over_rest_decay_width"]*width
            full=normed(S,H,[(16,6),(32,10)])
            f=full[-1]["rate_over_current_species_count_per_Hubble"]
            row=dict(T_GeV=d["T_GeV"],snapshot_sha256=r["snapshot_sha256"],full_support=full)
            if q==16:
                edges.append(S["lab_momentum_edges_GeV"][-1])
                assert min(f)>1e3
            else:
                assert max(f)<.003
                c=normed(crop(S,edges[j]),H,[(16,6),(32,10)])
                assert min(c[-1]["rate_over_current_species_count_per_Hubble"])>1e3
                row["same_q32_crop_both_species_to_q16"]=c
            result.append(row)
        rows.append(dict(q=q,stages=result))
    return dict(scope="Frozen open prompt comparison: heavy stationary decay and inverse only, "
                      "with empty exterior occupations but allowed final exterior states. "
                      "These one-way instantaneous rates are not sustained escape or actual "
                      "cohort formation, because depletion, re-entry, and feedback are absent.",
                cases=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tapes",type=Path,default=Path("results"))
    print(json.dumps(run(p.parse_args().tapes),indent=2,allow_nan=False))
