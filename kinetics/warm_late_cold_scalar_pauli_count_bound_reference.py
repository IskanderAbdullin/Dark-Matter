"""Bound open Pauli suppression by the full known light count, not a pointwise cap."""

import argparse
import json
from pathlib import Path

import numpy as np

from frw_exterior_restart_reference import import_restart
from three_species_decay_onsager_reference import M, P, m
from warm_late_cold_scalar_single_use_receipt_reference import run as reference
from warm_late_cold_scalar_single_use_receipt_reference import angular_vacuum_bound
from warm_light_occupation_stage_audit_reference import digest, read
from warm_prompt_reentry_shape_audit_reference import get_mom


def max_filled_length(N):
    pF=np.cbrt(6*np.pi**2*N/4)
    return pF**2/(np.hypot(m,pF)+m)


def test():
    p=np.r_[0.,np.geomspace(1e-8,6.,180)]
    Ea=np.hypot(m,p)
    n=4*np.diff(p**3)/(6*np.pi**2)
    # Use a difference that is still resolved for p close to zero.
    e=(p[1:]-p[:-1])*(p[1:]+p[:-1])/(Ea[1:]+Ea[:-1])
    rng=np.random.default_rng(92143)
    for slope in (0.2, .7, 5.):
        for _ in range(40):
            phase=rng.random(len(n))**slope
            assert phase@e<=max_filled_length(phase@n)*(1+6e-15)
    for j in (1,50,100,150,180):
        np.testing.assert_allclose(e[:j].sum(),max_filled_length(n[:j].sum()),rtol=3e-15)
    return dict(tested_random_positive_FD_profiles=120,exact_filled_step_controls=5)


def run(paths):
    proof=test()
    out=[]
    for path,z in zip(paths,reference(paths)["cases"],strict=True):
        x=read(path)
        S=x["snapshot"]
        assert z["hash"]==x["snapshot_sha256"]==digest(S)
        p=np.asarray(S["lab_momentum_edges_GeV"])
        n=get_mom(p,[S["entropy_variable_"+k] for k in ("scalar","fermion")])[0][1].sum()
        r=import_restart(S["positive_exterior_restart"]).array()
        out_count=r[5,r[1]==m].sum()
        limit=float(max_filled_length(n+out_count))
        rows=[]
        for a in z["controls"]:
            Hmin=a["E_parent_gate_GeV"]
            v=a["max_unpolarized_vacuum_angle_share"]
            assert abs(v-angular_vacuum_bound(a["q_over_T"]*S["T_visible_GeV"])[0])<1e-15
            fullwidth=2*P*np.sqrt((Hmin/M-1)*(Hmin/M+1))
            delta=limit/fullwidth
            assert 0<delta<1
            v_B=v*(1+a["archived_cold_BE_max_f"])/(1-delta)
            share=v_B*a["geometric_one_decay_upper_over_initial_cold_count"]
            rows.append(dict(p_over_T=a["q_over_T"],eligible_parent_minimum_open_width_GeV=float(fullwidth),
                relative_pauli_blocking_length_upper=float(delta),cold_decay_probability_upper=float(v_B),
                one_use_cold_supply_upper_over_initial_cold_count=float(share)))
        out.append(dict(archive=x["snapshot_sha256"],represented_light_number_GeV3=float(n),
                        tracked_light_number_GeV3=float(out_count),
                        maximum_occupied_energy_measure_GeV=limit,controls=rows))
    return dict(scope="Frozen archived scalar occupations in the cold band, open unpolarized decay and no "
                      "untracked light count. Pauli bound allows arbitrary isotropic light FD shape at known count. "
                      "The receipt population is used at most once. No claim about repeated use, polarization, "
                      "future Bose growth or an untracked light reservoir.",controls=proof,cases=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("archives",type=Path,nargs="+")
    args=p.parse_args()
    print(json.dumps(run(args.archives),indent=2,allow_nan=False))
