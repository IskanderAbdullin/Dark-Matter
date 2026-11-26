"""Audit the matched q32 and q64 trajectories with the nine-cell reverse lift."""

import argparse
import json
from pathlib import Path

import numpy as np

from frw_exterior_restart_reference import import_restart
from joint_prompt_two_moment_evolution_reference import Collision
from late_radiation_background_reference import MPL, TabulatedRadiation
from three_species_decay_onsager_reference import mp
from warm_coupled_serial_restart_reference import sha
from warm_history_lower_matched_audit_reference import measures
from warm_light_occupation_stage_audit_reference import digest, read
from warm_near_rest_light_inverse_decay_reference import window


def state(eos, spec, ref):
    p = np.asarray(spec["lab_momentum_edges_GeV"])
    a = spec["scale_factor"]
    a0, saved = ref["scale_factor"], ref["positive_exterior_restart"]
    b = a / a0
    ref_mom = measures(ref)
    moments = measures(spec)
    ex = import_restart(spec["positive_exterior_restart"])
    born = ex.array()
    inherited = import_restart(saved).array()
    k = inherited.shape[1]
    old = born[:,:k].copy()
    old[[0,5]] *= np.array([b,b**3])[:,None]
    np.testing.assert_allclose(old,inherited,rtol=1.6e-13,atol=0.)
    closure = abs((born[5,k:].sum()+moments[0].sum())*b**3/ref_mom[0].sum()-1.)
    assert closure < 8e-13
    T = spec["T_visible_GeV"]
    E = TabulatedRadiation(eos,T).state(1.)[1]+moments[1].sum()+ex.energy(1.)
    H = np.sqrt(E/3)/MPL
    operator = Collision(p)
    q = [np.asarray(spec["entropy_variable_"+key]) for key in ("scalar","fermion")]
    occ = moments[0,1,0]/operator.c[1].vol[0]
    prompt = operator.rate(q)[1,0,0]/(H*operator.c[1].vol[0])
    active = int(np.sum(operator.J==0))
    if np.hypot(mp,p[-1]) < window(p[1])[-2]:
        assert active==0 and prompt==0.
    return dict(T_GeV=T,relative_a=b,first_cell_occupation=float(occ),
                prompt_events=active,prompt_rate_per_H=float(prompt),count_defect=float(closure)), p,moments*a**3


def audit(eos, seed32, q32, wide, count):
    guide_tape = wide/"guide-profile.json"
    guide_rows = read(guide_tape)["blocks"]
    w0 = read(wide/"guide-start.json")
    r0 = read(seed32)[-1]
    manifest = read(wide/"correction-manifest.json")
    assert manifest["eos_sha256"] == sha(eos) and manifest["guide_sha256"] == sha(guide_tape)
    planned = len(guide_rows)
    assert planned in (24,48) and 0<=count<=planned
    initials = [x["snapshot"] for x in (r0,w0)]
    previous = [x["snapshot_sha256"] for x in (r0,w0)]
    D = measures(initials[0])*initials[0]["scale_factor"]**3
    scale = np.stack((D[0].sum(axis=1),np.full(2,D[1].sum()),np.full(2,D[1].sum())))
    out = []
    offset = None
    count_fallback = 0
    for i in range(0,count+1):
        if i:
            A = read(q32/f"step-{i}.json")
            B = read(wide/"q64"/f"step-{i}.json")
            assert A["index"]==B["step"]==i
            assert A==guide_rows[i-1]["checkpoint"]
            assert B["plan"]["guide_snapshots"][i-1]==A["snapshot_sha256"]
            assert B["plan"]["guides_sha256"]==manifest["guide_sha256"]
            assert A["previous"]==previous[0] and B["parent_snapshot_sha256"]==previous[1]
            assert A["plan"]["tape_sha256"]==sha(seed32)
        else:
            A,B=r0,w0
        row, arrays, rates = [],[],[]
        for j, v in enumerate((A,B)):
            S = v["snapshot"]
            previous[j]=v["snapshot_sha256"]
            assert digest(S)==previous[j]
            v0,p,u=state(eos,S,initials[j])
            assert abs(v0["relative_a"]/(1.10**(i/planned))-1.)<1.2e-13
            params=dict(snapshot_sha256=previous[j],**v0)
            if i:
                before=out[-1]["branches"][j]
                params["cold_slope_per_log_a"]=(params["first_cell_occupation"]-before["first_cell_occupation"])/np.log(params["relative_a"]/before["relative_a"])
                params["local_energy_defect"]=v["operator_summary"]["energy_ledger"]
                params["used_corrected_narrow_seed"]=v["audit"]["trial"]["used_q32_guess"] if j else None
                if j:
                    count_fallback += int(not params["used_corrected_narrow_seed"])
            row.append(params); arrays.append(p); rates.append(u)
        p,P=arrays;  x,y=rates
        assert abs(row[0]["T_GeV"]/row[1]["T_GeV"]-1.)<1.1e-13
        assert row[0]["prompt_events"]==0 and row[1]["prompt_events"]>0
        np.testing.assert_allclose(p,P[:len(p)],rtol=1.1e-13,atol=0.)
        signed = x-y[...,:x.shape[-1]]
        if offset is None:
            offset = signed
        L1 = abs(signed).sum(axis=-1)/scale
        drift = abs(signed-offset).sum(axis=-1)/scale
        frac=y[0,:,x.shape[-1]:].sum(axis=1)/y[0].sum(axis=1)
        out.append(dict(index=i,branches=row,shared_L1=L1.tolist(),
                        subsequent_difference_drift_L1=drift.tolist(),wide_only_count_fraction=frac.tolist()))
    lo, hi = out[0], out[-1]
    assert not any(x["branches"][0]["prompt_events"] for x in out)
    return dict(scope="Independent wider targets and scalar-pair support, narrow corrected states are initial guesses only.",
                fallback_steps=count_fallback, metric_order=["number","full_energy","pressure"],
                initial_shared_L1=lo["shared_L1"], endpoint_shared_L1=hi["shared_L1"],
                endpoint_subsequent_drift_L1=hi["subsequent_difference_drift_L1"],frames=out)


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("eos","seed32","q32","wide"):
        p.add_argument(key,type=Path)
    p.add_argument("count",type=int)
    a=p.parse_args()
    print(json.dumps(audit(a.eos,a.seed32,a.q32,a.wide,a.count),indent=2,allow_nan=False))
