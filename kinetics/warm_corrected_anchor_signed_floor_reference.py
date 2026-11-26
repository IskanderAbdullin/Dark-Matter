"""Check how matched corrected endpoints change an existing conditional floor."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import expit

from conditional_relic_normalization_reference import conversion
from late_radiation_background_reference import TabulatedRadiation
from three_species_decay_onsager_reference import m
from warm_light_occupation_stage_audit_reference import digest, read
from warm_prompt_reentry_shape_audit_reference import get_mom
from warm_signed_fermion_piecewise_clock_reference import assess


def run(eos,D,paths):
    first=read(D/"warm-isolated-fermion-survival-bound-2026-11-26.json")["cases"][-1]
    leptons=read(D/"warm-fermion-survival-with-leptons-2026-11-26.json")
    pions=read(D/"warm-signed-fermion-with-charged-pions-2026-11-26.json")["cases"][-1]["controls"][-1]
    early=TabulatedRadiation(eos,first["anchor_T_GeV"])
    k,Kpi=leptons["coefficient"],pions["pion_pair_sigma_v_cap_with_stated_support"]
    rows=[]
    for label,path in zip(("q32-12","q32-24","q32-48","q64-24"),paths,strict=True):
        row=read(path)
        S=row["snapshot"]
        assert digest(S)==row["snapshot_sha256"]
        T=S["T_visible_GeV"]
        rad=TabulatedRadiation(eos,T)
        b=np.exp((early.s0-rad.s0)/3)
        p=np.asarray(S["lab_momentum_edges_GeV"])
        theta=[np.asarray(S["entropy_variable_"+key]) for key in ("scalar","fermion")]
        n=get_mom(p,theta)[0][1].sum()
        O=m*n*conversion(eos)["Omega_h2_per_GeV_times_yield"]/(2*np.pi**2/45*np.exp(rad.s0))
        E=pions["new_available_energy_upper_GeV4"]/b**3
        data=assess(eos,T,n,E,k,Kpi,O,label)
        a,z=theta[1]
        maximum=float(np.max(expit(-np.minimum(a,a+z))))
        rows.append(dict(snapshot_sha256=row["snapshot_sha256"],
                         represented_light_maximum_at_anchor=maximum,
                         entropy_expansion_from_original_anchor=b,**data))
    np.testing.assert_allclose([r["T_GeV"] for r in rows],rows[0]["T_GeV"],rtol=2e-13)
    return dict(scope="Existing conditional selected-network comparison at matched numerical anchors. "
                      "Both fermion caps and the transported available-energy envelope remain hypotheses; "
                      "no untracked old number is credited and no extra energy is authorized.",
                branches=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    for x in ("q32_12","q32_24","q32_48","q64_24"):
        p.add_argument(x,type=Path)
    p.add_argument("--results",type=Path,default=Path("results"))
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.results,(a.q32_12,a.q32_24,a.q32_48,a.q64_24)),indent=2,allow_nan=False))
