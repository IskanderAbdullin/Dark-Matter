"""Audit sampled exterior receipts without subtracting whole elastic currents."""

import argparse
import json
from pathlib import Path

import numpy as np

from warm_light_occupation_stage_audit_reference import read


def run(paths):
    analyses=[read(path) for path in paths]
    assert all(len(t["rows"])==6 for t in analyses)
    assert len({tuple(t["order"]) for t in analyses})==len(analyses)
    out=[]
    for states in zip(*(t["rows"] for t in analyses),strict=True):
        a=states[0]
        assert all(x["source_hash"]==a["source_hash"] and x["time"]==a["time"] for x in states)
        assert all(x["physical_incoming_lab_ceiling_GeV"]==a["physical_incoming_lab_ceiling_GeV"] for x in states)
        ref=a["comparison_prompt_signed_over_rhoH"]
        assert abs(ref)>0.
        C=a["cut_GeV"]
        rows=[]
        for t in states:
            ext=np.asarray(t["newly_outside_by_species_N_E_W_over_rhoH"])
            n,e,w=ext.T
            np.testing.assert_allclose(e-C*n,w,rtol=3e-12,atol=0.)
            events,delta,Win=t["crossing_events_N_signed_W_incoming_W_over_rhoH"]
            assert all(t["old_and_new_shells"][j]>0 for j in (0,1))
            assert abs(sum(n)/events-1.)<3e-13
            assert sum(e)<= t["total_lab_energy_of_events_with_outgoing_exterior_over_rhoH"]*(1.+3e-13)
            assert abs(delta)<=sum(w)+Win+t["total_lab_energy_of_events_with_outgoing_exterior_over_rhoH"]
            own=np.asarray(t["old_chart_gross_and_signed"])
            nested=np.asarray(t["matched_chart_closed_gross_signed"])
            assert max(abs(own-nested))/max(own[:2])<2e-9
            rows.append(dict(newly_outgoing_W=sum(w),
                     crossing_event_W_change=delta,relative_to_absolute_prompt=max(sum(w),abs(delta))/abs(ref),
                     closed_signed=own[2],newly_outgoing_energy=sum(e)))
        out.append(dict(hard_GeV=a["hard_GeV"],time=a["time"],input_hash=a["source_hash"],prompt=ref,controls=rows))
    return dict(scope="Pairwise ledger audits; these maxima are over selected states and quadrature, not rigorous upper bounds.",
                per_state=out,maximum_normalized_exterior_W=max(r["newly_outgoing_W"] for t in out for r in t["controls"]),
                max_prompt_relative=max(r["relative_to_absolute_prompt"] for t in out for r in t["controls"]))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("tapes",type=Path,nargs="+")
    a=p.parse_args()
    print(json.dumps(run(a.tapes),indent=2,allow_nan=False))
