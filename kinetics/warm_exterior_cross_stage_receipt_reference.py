"""Audit one-way exterior receipts across a serial plan boundary."""
import argparse
import json
from pathlib import Path

from warm_positive_exterior_receipt_reference import read,run_case
from warm_observed_suboperator_audit_reference import audit
from warm_light_occupation_stage_audit_reference import verify


def run(pions,peripheral,exterior,parents,earlier,later):
    P,B,X=map(read,(pions,peripheral,exterior))
    result=[]
    for a,b,x,seed,old,new in zip(P["cases"],B["cases"],X["cases"],parents,earlier,later,strict=True):
        q=a["q"]
        assert q==b["q"]==x["q"]
        verify(seed,old)
        audit(new,old)
        Q=read(old)+[z["checkpoint"] for z in read(new)["blocks"]]
        assert len(Q)==6
        E=a["controls"][-1]["new_available_energy_upper_GeV4"]
        R=b["max_gross_comoving_peripheral_input_fraction_for_test"]
        rows=run_case(seed,Q,E,R,x["old_exterior_rho_pressure_only_upper_GeV4"],
                      require_sequential_labels=False)
        # This calculation never imposes a source outside the one-way
        # model; its output is a cumulative receipt, not reacting use.
        oldnew=rows[-1]["cumulative_high_gross_created_GeV4"]/rows[3]["cumulative_high_gross_created_GeV4"]-1.
        result.append(dict(q=q,relative_growth_of_receipt_after_later_plan=oldnew,rows=rows))
    return dict(scope="Positive high exterior creation continued across two serial plans. "
                      "No exterior reactions modeled. Older energy upper is conditional; "
                      "not an inference of cumulative peripheral reaction use.",
                cases=result)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("pions","peripheral","exterior"):
        p.add_argument(k)
    for k in ("parents","earlier","later"):
        p.add_argument("--"+k,nargs=2,required=True,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.pions,a.peripheral,a.exterior,a.parents,a.earlier,a.later),
                     indent=2,allow_nan=False))
