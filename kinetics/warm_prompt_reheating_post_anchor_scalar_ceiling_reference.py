"""Weighted-count ceiling allowing selected dark heavy regeneration after cold anchors."""
import argparse
import json
from pathlib import Path
import numpy as np

from archived_collisionless_streaming_reference import Background
from conditional_relic_normalization_reference import conversion
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import M,mp
from warm_prompt_reheating_expanding_slow_reference import visible


def run(eos,history,kernels):
    a=json.loads(history.read_text())
    k=json.loads(kernels.read_text())
    C=k["simultaneous_energy_density_coefficient_GeVminus2"]
    b=Background(eos)
    unit=conversion(eos)["Omega_h2_per_GeV_times_yield"]*mp
    out=[]
    for row in (a["arms"][-1]["trajectory"][0],a["arms"][-1]["trajectory"][-1]):
        tv=row["T_visible"]
        _,s=visible(b,tv)
        ed=sum(row["energy_per_visible_entropy"])*s
        H=np.sqrt(np.pi**2/45)*tv*tv/MPL
        q=sum(row["yields"])
        delta=C*ed*ed/(M*M*H*s)
        out.append(dict(visible_T=tv,initial_weighted_yield=q,
                        initial_dark_energy_density=ed,photon_H_floor_at_anchor=H,
                        extra_weighted_yield_upper_if_isolated=delta,
                        scalar_Omega_h2_upper_if_all_counts_are_converted=unit*(q+delta),
                        scalar_share_of_reference_half=unit*(q+delta)/.060))
    return dict(scope="Isolated future with no added energy or particles. Allows arbitrary elastic redistribution, prompt "
                      "decay and its reverse, selected light and scalar pair reactions and heavy regeneration, using "
                      "the independently bounded vacuum kernel. No other weighted count source; no condensate "
                      "dynamics, external heating or new interactions. Photon expansion floor and no exotic "
                      "visible entropy. The chosen input histories remain controlled diagnostics, not a solved relic.",
                bounded_heavy_coefficient=C,anchors=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for x in ("eos","history","kernels"):p.add_argument(x,type=Path)
    a=p.parse_args();print(json.dumps(run(a.eos,a.history,a.kernels),indent=2,allow_nan=False))
