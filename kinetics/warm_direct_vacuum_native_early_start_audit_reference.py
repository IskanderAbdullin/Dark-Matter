"""Audit the controlled intervention on CLASS's opacity-threshold start."""
import argparse,json
from pathlib import Path
import numpy as np
from warm_direct_vacuum_native_focus_audit_reference import run as audit


def run(old,middle,early,pilot):
    phases=[audit(p,pilot) for p in (old,middle,early)]
    response=np.array([[i["exact_k_physical_minus_trace"] for i in a["cases"]] for a in phases])
    offset=np.array([[i["exact_k_trace_presence_offset"] for i in a["cases"]] for a in phases])
    early_starts={a["case"]:np.asarray(a["first_tau_Mpc_and_scale_at_exact_modes"]) for a in phases[-1]["initial_times"]}
    assert np.max(abs(early_starts["trace"][:,1]/early_starts["base"][:,1]-1))<1e-9
    assert np.max(abs(early_starts["physical"][:,1]/early_starts["trace"][:,1]-1))<1e-7
    assert np.max(abs(offset[2]))<2e-8 and np.max(abs(offset[0]))>1e-5
    assert np.max(abs(response[2]-response[0]))<2e-10
    return dict(scope="Exact-mode RK comparison with no change to species PSD, density, or hierarchy. Only "
                      "the opacity initialization threshold changes. Evidence for the start-time contribution "
                      "to the zero-species artifact; no claim about all CLASS approximation branches.",
                opacity_thresholds=["CLASS default",3e-6,2e-7],
                max_zero_vs_trace_offsets=np.max(abs(offset),axis=(1,2)).tolist(),
                maximum_tight_start_shift_in_physical_minus_trace=float(np.max(abs(response[2]-response[0]))),
                results=phases)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("default","middle","early","pilot"):
        p.add_argument(k,type=Path)
    a=p.parse_args();print(json.dumps(run(a.default,a.middle,a.early,a.pilot),indent=2,allow_nan=False))
