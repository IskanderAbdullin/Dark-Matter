"""Compare portable metric integration to the archived raw-metric results."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run(root):
    keep=[]
    for s in ("off","on"):
        paths=[root/f"finite-cooling-los-{s}-{i}-2026-11-26.json" for i in ("raw-reference","portable","fixture")]
        raw,port,manifest=[json.loads(p.read_text()) for p in paths]
        assert raw["eos_sha256"]==port["eos_sha256"] and raw["spectrum_sha256"]==port["spectrum_sha256"]
        for a,b,c in zip(raw["results"],port["results"],manifest["files"]):
            assert a["k_Mpc_minus1"]==b["k_Mpc_minus1"]==c["k_Mpc_minus1"]
            metric=root/f"los-metric-{s}"/c["file"]
            assert a["metric_sha256"]==b["original_native_metric_sha256"]==c["raw_metric_sha256"]
            assert b["metric_sha256"]==c["fixture_sha256"]==sha(metric)
            r,t=[v["cases"][-1] for v in (a,b)]
            key="los_over_cdm_joined_prompt_recoil_losIBP_direct_monopoleIC"
            rho=float(np.max(abs(np.asarray(r[key])-np.asarray(t[key]))))
            key2="independent_deltaP_over_deltaRho_joined_prompt_recoil"
            pres=float(np.max(abs(np.asarray(r[key2])/np.asarray(t[key2])-1)))
            assert rho<3e-8 and pres<2e-6
            keep.append(dict(metric_driver=s,k_h_per_Mpc=b["k_h_per_Mpc"],fixture_bytes=metric.stat().st_size,
                 max_abs_density_compaction_change=rho,max_relative_pressure_compaction_change=pres))
    return dict(rows=keep,scope="Comparison of archived raw-driver results and compressed-driver outputs; the source metric is prescribed, not reconstructed or evolved.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("results",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.results),indent=2,allow_nan=False))
