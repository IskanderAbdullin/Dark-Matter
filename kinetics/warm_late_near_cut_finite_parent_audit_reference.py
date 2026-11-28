"""Separate finite-parent effects from taking fewer implicit substeps."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from warm_light_occupation_stage_audit_reference import digest, read


def run(matched, finite, fine):
    cases=[]
    for a,b,c in zip(matched,finite,fine,strict=True):
        U=read(a/"summary.json")
        V=read(b/"summary.json")
        R=read(c/"summary.json")
        t=U["rows"]
        f=V["rows"]
        r=R["rows"]
        u=U["plan"]
        assert u["seed_sha256"]==V["initial_hash"]==R["plan"]["seed_sha256"]
        assert V["initial_parent"]=="equilibrium" and u["substeps"]==1 and R["plan"]["substeps"]==3
        assert [j["time"] for j in t[1:]]==[j["time"] for j in f]==[j["time"] for j in r[1:]]
        assert V["eos_sha256"]==u["eos_sha256"]
        w=float(t[0]["tail_energy"])
        records=[]
        for i,(x,y,z) in enumerate(zip(t[1:],f,r[1:],strict=True),start=1):
            assert y["seeding_phase_hash"]==z["hash"] and len(x["substep_solver"])==1
            s=read(b/f"step-{i}.json")
            h=np.asarray(s["parent_occupations"],dtype="<f8")
            assert len(h)==V["parent_energy_nodes"] and np.all((h>=0)&(h<=1))
            assert s["hash"]==digest(s["snapshot"]) and s["row"]==y
            assert max(abs(np.asarray(y["total_invariant_relative_defects"])))<2e-12
            assert abs(y["stationary_W"]/z["tail_energy"]-1.)<4e-11
            assert abs(y["actual_W"]-x["tail_energy"])<1e-3*w
            rates=[float(v) for v in (x["signed_native_flux"],y["actual_signed_flux"],z["signed_native_flux"])]
            assert all(np.sign(rates[0])==np.sign(q) for q in rates)
            records.append(dict(time=x["time"],W_stationary_one=x["tail_energy"],W_finite_parent=y["actual_W"],
                W_stationary_three=z["tail_energy"],matched_parent_shift_over_initial_W=float((y["actual_W"]-x["tail_energy"])/w),
                rates_stationary_one_finite_stationary_three=rates,parent_phase_sha256=hashlib.sha256(h.tobytes()).hexdigest(),
                heavy_energy_change_over_rho=y["signed_heavy_energy_change_over_rho"], daughter_sha256=s["hash"],
                total_invariant_relative_defects=y["total_invariant_relative_defects"],
                parent_equation_roundoff_units=y["parent_weighted_residual_over_binary64_rounding_scale"]))
        cases.append(dict(halfwidth_GeV=u["shell_halfwidth_GeV"],initial_W=w,source_sha256=u["seed_sha256"],
                          parent_nodes=V["parent_energy_nodes"], rows=records))
    return dict(scope="Fixed background, closed prompt. Matched one-substep comparison separates parent lifetime "
                "from the stationary one-versus-three numerical clock change.",cases=cases)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("matched","finite","fine"):
        p.add_argument("--"+k,nargs=2,type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.matched,a.finite,a.fine),indent=2,allow_nan=False))
