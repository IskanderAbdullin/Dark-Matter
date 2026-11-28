"""Time-ladder audit after the conservative nine-cell lift on q32."""

import argparse
import json
from pathlib import Path

import numpy as np

from frw_exterior_restart_reference import import_restart
from warm_coupled_serial_restart_reference import sha
from warm_history_lower_matched_audit_reference import measures
from warm_light_occupation_stage_audit_reference import digest, read


def series(eos, seed, folder, planned, limit):
    context = read(folder/"cold_partition_manifest.json")
    assert context["cut_cells"]==9 and context["steps"]==planned
    a0w = read(seed)[-1]
    A = a0w["snapshot"]; a0 = A["scale_factor"]
    last = a0w["snapshot_sha256"]
    assert digest(A)==last
    e0 = measures(A)
    old = import_restart(A["positive_exterior_restart"]).array()
    temperature, moment, occupation, stats = [], [], [], []
    for i in range(limit+1):
        if i:
            row = read(folder/f"step-{i}.json")
            assert row["previous"] == last and row["index"] == i
            assert row["plan"]["tape_sha256"]==sha(seed) and row["plan"]["EOS_sha256"]==sha(eos)
            A = row["snapshot"]
            last = row["snapshot_sha256"]
            assert digest(A) == last
        a = A["scale_factor"]
        b = a/a0
        assert abs(b/(1.1**(i/planned))-1.)<1e-13
        p = np.asarray(A["lab_momentum_edges_GeV"])
        m = measures(A)
        ext = import_restart(A["positive_exterior_restart"]).array()
        n = old.shape[1]
        z = ext[:,:n].copy()
        z[[0,5]] *= np.asarray([b,b**3])[:,None]
        np.testing.assert_allclose(z,old,rtol=1.6e-13,atol=0.)
        num = abs((m[0].sum()+ext[5,n:].sum())*b**3/e0[0].sum()-1.)
        assert num<8e-13
        temperature.append(A["T_visible_GeV"])
        moment.append(m*a**3)
        occupation.append(float(m[0,1,0]/(4*p[1]**3/(6*np.pi**2))))
        if i:
            stats.append(dict(step=i,snapshot_sha256=last,number_defect=float(num),
                          local_energy_defect=row["operator_summary"]["energy_ledger"],
                          own_history=row["audit"]["solver"]["used_history"]))
    slope = np.diff(occupation)/(np.log(1.10)/planned)
    return dict(planned_steps=planned,count=limit,T_GeV=temperature,first_cold_occupation=occupation,
                cold_slope_per_log_a=slope.tolist(),stats=stats),np.asarray(moment)


def audit(eos,seed,refs,fin,planned=24,limit=12):
    r, a = series(eos,seed,refs,planned//2,limit)
    s, b = series(eos,seed,fin,planned,2*limit)
    np.testing.assert_allclose(r["T_GeV"],s["T_GeV"][::2],rtol=1e-13)
    norm = a[0].sum(axis=-1)
    L1 = abs(a-b[::2]).sum(axis=3)/norm
    ra = np.asarray(r["cold_slope_per_log_a"])
    rb = np.asarray(s["cold_slope_per_log_a"]).reshape(limit,2).mean(axis=1)
    return dict(scope="Same q32 inputs and corrected operator, only expansion time grid refined. One comparison does not prove an asymptotic order.",
                metric_order=["number","energy","pressure"],branches=[r,s],trajectory_L1=L1.tolist(),
                endpoint_L1=L1[-1].tolist(),same_interval_first_cold_slope_difference=(ra-rb).tolist())


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","seed","coarse","fine"):
        p.add_argument(k,type=Path)
    p.add_argument("--factor-count",type=int,default=24)
    p.add_argument("--limit",type=int,default=12)
    a=p.parse_args()
    print(json.dumps(audit(a.eos,a.seed,a.coarse,a.fine,a.factor_count,a.limit),
                     indent=2,allow_nan=False))
