"""Analyze pair rates on the extended warm state and its exact inner crop."""
import argparse
import json
from pathlib import Path

import numpy as np

from archived_scalar_pair_reaction_reference import run
from sm_dark_entropy_exchange_reference import StandardModel


def digest(result):
    a=result["archives"]
    assert len(a)==3
    assert "same-state crop" in a[2]["snapshot"]
    assert max(abs(x["T_GeV"]/a[0]["T_GeV"]-1) for x in a)<2e-12
    rows=[]
    for v in a:
        r=v["rows"]
        k=r[-1]
        F=k["conditional_net_stable_number_loss_per_Hubble"]/(1-k["reverse_over_stimulated_forward"])
        assert F>0. and k["relative_population_quadrature_defect"]<.02
        rows.append(dict(snapshot=v["snapshot"],T_GeV=v["T_GeV"],order=k["order"],
            normalized_gross_stable_loss_from_forward=F,
            normalized_gross_stable_gain_from_reverse=F*k["reverse_over_stimulated_forward"],
            normalized_signed_net_stable_loss=k["conditional_net_stable_number_loss_per_Hubble"],
            stable_inverse_over_stimulated_forward=k["reverse_over_stimulated_forward"],
            scalar_retained_number_and_energy_fraction=[k["retained_fraction_of_produced_scalar_number"],
                                                          k["retained_fraction_of_produced_scalar_energy"]],
            stimulated_over_vacuum_forward=k["stimulated_over_vacuum_forward"],
            coarse_to_fine_net_loss_relative_change=float(
                k["conditional_net_stable_number_loss_per_Hubble"]/r[0]["conditional_net_stable_number_loss_per_Hubble"]-1.),
            normalized_pair_ledger_defects=[v["pair_number_ledger_defect"],v["opposite_energy_ledger_defect"]]))
    keys=("normalized_gross_stable_loss_from_forward","normalized_gross_stable_gain_from_reverse",
          "normalized_signed_net_stable_loss")
    return dict(rows=rows,quantities_compared=keys,
          broad_over_same_state_crop_minus_one=[rows[1][k]/rows[2][k]-1. for k in keys],
          same_state_crop_over_original_run_minus_one=[rows[2][k]/rows[0][k]-1. for k in keys],
          scope="Signed tree pair operator on archived fixed 12 MeV states only. Radiation H, normalization per retained stable count of each state, exact inner crop. No time evolution or proof of complete outer convergence.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos"); p.add_argument("original_shape"); p.add_argument("wide_shape")
    p.add_argument("--reuse-native-rates",type=Path)
    a=p.parse_args()
    if a.reuse_native_rates:
        raw=json.loads(a.reuse_native_rates.read_text())
    else:
        e=StandardModel(a.eos,lower_knot_limit=.009)
        raw=dict(archives=run(e,[a.original_shape,a.wide_shape],[(3,7,6,2),(5,10,8,3)],16.))
    print(json.dumps(dict(summary=digest(raw),raw=raw),indent=2,allow_nan=False))
