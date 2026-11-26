"""Cross-check inclusive low-reheat source and scalar-pair histories."""
import argparse
import json
from pathlib import Path

import numpy as np

from conditional_relic_normalization_reference import conversion, proxy
from sm_dark_entropy_exchange_reference import StandardModel
from visible_frw_energy_ledger_reference import ConservativeEOS


def rows(path):
    return json.loads(Path(path).read_text())["rows"]


def counts(row,eos,factor):
    s0=eos.state(row["initial_visible_T_GeV"])[2]
    v=row["visible_bath_feedback"]
    sf=s0+v["visible_comoving_entropy_change"]
    kept=np.array(row["final_counts_div_initial_visible_entropy"])*s0/sf
    outside=(np.array(row["comoving_source_channel_count_and_total_energy_tails"])
            +np.array(row["first_pulse_count_and_total_energy_tails"]))[:,0].sum(axis=1)/sf
    tail=np.array([outside[2],outside[:2].sum()])
    delta=np.asarray(row.get("pair_delta_retained_counts_relative_to_final_retained",[0.,0.]))*kept
    pair_tail=np.asarray(row.get("pair_exterior_counts_relative_to_final_retained",[0.,0.]))*kept
    scale=max(np.max(abs(delta)),1e-280)
    assert abs(sum(delta+pair_tail))/scale<2e-9
    result=kept+tail+pair_tail
    return dict(steps=row["nsteps"],enabled=bool(row.get("included_scalar_pair_conversion",False)),
        signed_visible_portal=row["state_dependent_quantum_visible_source"],
        final_visible_temperature_GeV=row["final_visible_temperature_GeV"],
        retained_fixed_Omega_h2=proxy(kept,factor),
        inclusive_fixed_Omega_h2=proxy(result,factor),
        exterior_fixed_Omega_h2=proxy(tail+pair_tail,factor),
        pair_production_and_loss_yield=(delta+pair_tail).tolist(),
        integrated_pair_contribution_relative_to_final_retained=(delta/kept).tolist(),
        pair_peak_inverse_over_forward=row.get("pair_max_inverse_over_forward",0.),
        reported_number_defect=row["total_number_defect"],
        visible_dark_exterior_energy_defect=v.get("visible_dark_exterior_joint_ledger_defect_in_initial_visible_units"))


def comparison(table,inputs,late,bridge_input_index=1):
    eos=ConservativeEOS(StandardModel(table,lower_knot_limit=.001))
    fac=conversion(table)["Omega_h2_per_GeV_times_yield"]
    chains=[]
    for filename in inputs:
        rr=rows(filename)
        q=[counts(r,eos,fac) for r in rr]
        chains.append(dict(input=Path(filename).name,rows=q))
    if not late:
        return chains,None
    q=chains[bridge_input_index]["rows"][-1]
    assert q["enabled"]
    t0=q["final_visible_temperature_GeV"]
    r=json.loads(Path(late).read_text())["states"][-1]
    assert abs(r["initial_and_final_visible_temperatures_GeV"][0]/t0-1.)<1e-9
    a=np.array(q["retained_fixed_Omega_h2"])
    b=np.array(q["exterior_fixed_Omega_h2"])
    c=a*np.array(r["inclusive_species_counts_over_initial"])
    return chains,dict(early_fixed_count_Omega_h2=(a+b).tolist(),
        late_updated_fixed_count_Omega_h2_with_early_exterior_carried=(c+b).tolist(),
        late_change_relative_to_early_combined=((c+b)/(a+b)-1.).tolist(),
        late_scale_factor=r["scale_factor"],
        note="Late stage evolves only represented scalar-pair reactions on prescribed radiation. Early exterior counts carried collisionlessly without re-entry.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("table")
    p.add_argument("inputs",nargs=3,help="matched/source, newly paired vacuum, and newly paired quantum")
    p.add_argument("--late")
    p.add_argument("--bridge-input-index",type=int,choices=[1,2],default=1)
    a=p.parse_args()
    x,z=comparison(a.table,a.inputs,a.late,a.bridge_input_index)
    print(json.dumps(dict(controls=x,bridge=z),indent=2,allow_nan=False))
