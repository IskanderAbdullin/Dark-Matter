"""Re-evaluate partial pointwise direct light gain on the late corrected clock."""

import argparse
import json
from pathlib import Path

import numpy as np

from warm_heavy_band_aggregate_energy_bound_reference import inverse_bound
from three_species_decay_onsager_reference import M, P
from warm_isolated_recycling_energy_envelope_reference import radiation_bound
from warm_light_direct_pair_dose_reference import direct_dose
from warm_light_occupation_stage_audit_reference import digest,read
from warm_signed_fermion_piecewise_clock_reference import weights


def run(eos,D,anchor):
    cap=read(D/"warm-signed-fermion-with-charged-pions-2026-11-26.json")["cases"][-1]["controls"][-1]
    endpoint=read(anchor)
    assert endpoint["snapshot_sha256"]==digest(endpoint["snapshot"])
    rows=read(D/"warm-corrected-anchor-signed-floor-2026-11-26.json")["branches"]
    row=next(r for r in rows if r["snapshot_sha256"]==endpoint["snapshot_sha256"])
    T=endpoint["snapshot"]["T_visible_GeV"]
    assert row["T_GeV"]==T
    E=row["energy_upper_GeV4"]
    Kpi=cap["pion_pair_sigma_v_cap_with_stated_support"]
    f0=row["represented_light_maximum_at_anchor"]
    # The available slack here is for a represented initial mode.
    # This does not claim that unrecorded modes start below f0.
    W=radiation_bound(eos,T,500.)
    b,h,H=W["b_stop"],W["floor_min"],W["photon_anchor_GeV"]
    powers=(3,4,6)
    uniform={p:((1-b**(2-p))/h+b**(2-p))/(H*(p-2)) for p in powers}
    loc=weights(eos,T,b,2800,H,powers=powers).sum(axis=0)
    ratios=loc/np.array(list(uniform.values()))
    assert np.max(ratios)<1+3e-10
    local=dict(zip(powers,loc,strict=True))
    coarse=direct_dose(eos,T,E,Kpi,uniform.__getitem__,margin=.5-f0)
    fine=direct_dose(eos,T,E,Kpi,local.__getitem__,margin=.5-f0)
    frac=fine["diagnostic_heavy_energy_fraction_for_direct_dose_half"]
    assert frac>=coarse["diagnostic_heavy_energy_fraction_for_direct_dose_half"]
    for d in (coarse,fine):
        d["residual_after_represented_initial_occupation_and_direct_dose"]=d.pop(
            "residual_margin_to_half_ignoring_initial_and_other_sources")
    old_b=row["entropy_expansion_from_original_anchor"]
    old=read(D/"warm-heavy-band-aggregate-energy-bound-2026-11-26.json")["cases"][-1]
    bound=next(r for r in old["rows"] if r["light_cap"]==.5)
    barrier=read(D/"warm-heavy-internal-source-barrier-2026-11-26.json")["results"][-1]
    assert old["q"]==barrier["q"]
    extra=.25*next(r for r in barrier["cases"] if r["cap"]==.5)[
        "total_bounded_input_over_outgoing_decay_margin"]
    bands=[]
    for c in bound["controls"]:
        oldcut=c["cut_GeV"]
        cut=oldcut/old_b
        initial=c["strongest_bound_for_below_cut"]
        def future(n,k):
            u,*_=inverse_bound(E,.5,cut,n,k)
            return float(2*np.hypot(M,cut)/E*(M*u/(np.pi**2*P)+2*cut**3*extra/(3*np.pi**2)))
        local_band=[future(n,k) for n,k in ((400,40),(1200,120),(2400,240))]
        X=max(initial,local_band[-1])
        bands.append(dict(cut_at_original_anchor_GeV=oldcut,cut_here_GeV=cut,
                          older_inherited_cap_if_assumptions_hold=initial,
                          local_source_over_removal_cap_only=local_band,
                          below_cut_cap_including_inherited_cap=X,
                          remaining_above_cut_fraction_for_direct_subtotal=frac-X))
    return dict(scope="Conditional partial direct pair gain only. Assumes heavy cap 1/2 and the "
                      "stated transported gross-energy budget. No light decay/inverse or elastic gain "
                      "is allowed in this subtotal; no cap for the omitted modes or future solution follows.",
                snapshot_sha256=row["snapshot_sha256"], T_GeV=T, energy_upper_GeV4=E,
                represented_light_maximum_at_anchor=f0, represented_mode_slack_to_half=.5-f0,
                powers=powers,local_weight_over_uniform_floor=ratios.tolist(),
                uniform_clock=coarse,piecewise_clock=fine,
                transported_band_scope="Initial below-cut cap is admitted only under all hypotheses "
                "of the original bound, including initially empty old heavy band and the earlier light "
                "cap. Uses original conservative extra source; no new unknown transfer into the band. "
                "Above-cut allowance is a target for the partial direct subtotal, not an upper bound.",
                transported_finite_bands=bands)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("anchor",type=Path)
    p.add_argument("--results",type=Path,default=Path("results"))
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.results,a.anchor),indent=2,allow_nan=False))
