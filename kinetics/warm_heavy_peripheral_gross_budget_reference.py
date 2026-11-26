"""Translate the above-cut heavy allowance into gross incoming peripheral energy."""
import argparse
import json
from pathlib import Path

from three_species_decay_onsager_reference import m,M,mp


def run(high_path):
    old=json.loads(Path(high_path).read_bytes())
    out=[]
    for a in old["cases"]:
        pc=a["high_cut_GeV"]; p=a["assumed_equal_comoving_daughter_cap_float_GeV"]
        assert max(a["internal_produced_heavy_cap_GeV"])<pc
        ratio=max(p/pc,m/M,mp/M)
        slack=a["light_direct_diagnostic_limit"]-a["combined_heavy_energy_fraction_upper"]
        assert 0<slack<1 and 0<ratio<1
        out.append(dict(q=a["q"],comoving_high_cut_GeV=pc,
                        core_daughter_energy_over_high_parent_max=ratio,
                        needed_min_outside_energy_share_per_event=1-ratio,
                        remaining_high_H_fraction_budget=slack,
                        max_gross_comoving_peripheral_input_fraction_for_test=slack*(1-ratio)))
    return dict(scope="Conditional budget on the gross received energy from above-cut "
                      "light/scalar inputs to high heavy production. Does not bound actual "
                      "gross usage from an instantaneous or net exterior reservoir, nor "
                      "include elastic heavy transfer across the cut.",cases=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("high")
    a=p.parse_args()
    print(json.dumps(run(a.high),indent=2,allow_nan=False))
