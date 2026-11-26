"""Compare matched warm source histories after widening only the outer bins."""
import argparse
import json
from pathlib import Path

import numpy as np

from conditional_relic_normalization_reference import conversion
from early_scalar_pair_budget_reference import counts
from sm_dark_entropy_exchange_reference import StandardModel
from visible_frw_energy_ledger_reference import ConservativeEOS
from scalar_pair_expanding_continuation_reference import moments,MASS
from scipy.special import expit


def last(path,kind):
    return json.loads(Path(path).read_text())[kind][-1]


def run(eosfile,standard,standard_shape,wide,wide_shape):
    eos=ConservativeEOS(StandardModel(eosfile,lower_knot_limit=.001))
    f=conversion(eosfile)["Omega_h2_per_GeV_times_yield"]
    z=[last(x,"rows") for x in (standard,wide)]
    S=[last(x,"snapshots") for x in (standard_shape,wide_shape)]
    support=[]
    for r,s in zip(z,S):
        assert r["nsteps"]==s["steps"] and r["nsteps"]==48
        assert r["included_scalar_pair_conversion"] and r["state_dependent_quantum_visible_source"]
        assert r["below_pole_elastic"] and s["below_pole_elastic"]
        assert abs(r["final_a"]/s["scale_factor"]-1)<2e-12
        assert abs(r["final_visible_temperature_GeV"]/s["T_visible_GeV"]-1)<2e-12
        q=np.asarray(s["lab_momentum_edges_GeV"])*s["scale_factor"]/s["T_initial_GeV"]
        support.append(q)
    a,b=support
    assert len(b)==len(a)+5 and np.max(abs(a-b[:len(a)]))<3e-12
    rs=[counts(r,eos,f) for r in z]
    outer_start=len(a)-1
    th=[np.asarray(S[1][k]) for k in ("entropy_variable_scalar","entropy_variable_fermion")]
    mom=moments(np.asarray(S[1]["lab_momentum_edges_GeV"]),th)
    energy=mom[:,1]+MASS[:,None]*mom[:,0]
    cen=[float(np.min(np.minimum(x[0,outer_start:],x.sum(axis=0)[outer_start:]))) for x in th]
    assert cen[0]>0
    tail_population=dict(number_fraction_by_species=(mom[:,0,outer_start:].sum(axis=1)/mom[:,0].sum(axis=1)).tolist(),
                total_energy_fraction_by_species=(energy[:,outer_start:].sum(axis=1)/energy.sum(axis=1)).tolist(),
                exact_outer_cell_max_occupations=[float(np.exp(-cen[0])/(-np.expm1(-cen[0]))),float(expit(-cen[1]))])
    def vec(k):
        return [q[k] for q in rs]
    key="inclusive_fixed_Omega_h2"
    A,B=vec(key)
    em=[r["visible_bath_feedback"]["exterior_a4_energy_one_mean_per_cohort"] for r in z]
    en=[r["visible_bath_feedback"]["visible_dark_exterior_joint_ledger_defect_in_initial_visible_units"] for r in z]
    return dict(order=["original q/T_R up to 16","same inner bins plus five to 32"],
                cases=rs,relative_wide_vs_original_inclusive_Omega=(np.asarray(B)/A-1).tolist(),
                number_of_original_and_extended_bins=[len(a)-1,len(b)-1],
                matched_inner_edges_max_absolute_difference=float(np.max(abs(a-b[:len(a)]))),
                outside_a4_energy=em,outside_energy_wide_over_original=em[1]/em[0],
                joint_energy_ledger_defects=en, new_resolved_outer_population=tail_population,
                scope="Matched full signed visible-source, prompt, elastic and pair evolution from the initially empty 50 MeV benchmark to 12 MeV. Extra outer support at fixed inner and time resolution. Fixed-count proxies, not final abundances or all-support convergence.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for n in ("eos","original","original_shape","wide","wide_shape"):
        p.add_argument(n)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.original,a.original_shape,a.wide,a.wide_shape),indent=2,allow_nan=False))
