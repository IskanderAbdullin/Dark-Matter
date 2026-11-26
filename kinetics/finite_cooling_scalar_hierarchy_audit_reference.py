"""Cross-check identical scalar input under two native hierarchy settings."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

from archived_class_scalar_tracer_reference import density_ratio


def fields(path):
    return {k.strip():v.strip() for k,v in (l.split("=",1) for l in path.read_text().splitlines()
                                          if l.strip() and not l.lstrip().startswith("#"))}


def run(root,tags):
    rows=[]
    for tag in tags:
        normal=f"finite_cooling_tracer_{tag}_s10_sync"
        strict=f"finite_cooling_tracer_{tag}_s10_hierarchy_sync"
        a,b=[fields(root/(f+".ini")) for f in (normal,strict)]
        assert a.keys()==b.keys()
        changed={k for k in a if a[k]!=b[k]}
        assert changed=={"ncdm_fluid_trigger_tau_over_tau_k","l_max_ncdm","root"}
        assert (int(a["l_max_ncdm"]),int(b["l_max_ncdm"]))==(14,22)
        x=[]
        for i,z in enumerate((0,5,20),1):
            q0,p0=density_ratio(root/f"{normal}_00_z{i}_tk.dat",z)
            q1,p1=density_ratio(root/f"{strict}_00_z{i}_tk.dat",z)
            assert min(q0[-1],q1[-1])>=10. and max(q0[0],q1[0])<=.015
            grid=np.unique(np.r_[np.geomspace(.015,10.,1001),q0[(q0>.015)&(q0<10.)],q1[(q1>.015)&(q1<10.)]])
            f=[PchipInterpolator(np.log(q),p) for q,p in ((q0,p0),(q1,p1))]
            dp=f[1](np.log(grid))-f[0](np.log(grid))
            j=np.argmax(np.abs(dp))
            probes=np.array([.015,.03,.1,.3,1.,2.,4.,8.,10.])
            x.append(dict(z=z,largest_absolute_increment=float(abs(dp[j])),
                    signed_increment_at_largest_probe=float(dp[j]),k_at_largest_increment_h_per_Mpc=float(grid[j]),
                    k_probes_h_per_Mpc=probes.tolist(),
                    baseline_delta_scalar_over_delta_cdm=f[0](np.log(probes)).tolist(),
                    stricter_delta_scalar_over_delta_cdm=f[1](np.log(probes)).tolist()))
        rows.append(dict(tracer_history=tag,baseline=normal,stricter=strict,
             PSD_sha256=hashlib.sha256(Path(a["ncdm_psd_filenames"]).read_bytes()).hexdigest(),
             baseline_hierarchy_and_trigger=[int(a[k]) for k in ("l_max_ncdm","ncdm_fluid_trigger_tau_over_tau_k")],
             stricter_hierarchy_and_trigger=[int(b[k]) for k in ("l_max_ncdm","ncdm_fluid_trigger_tau_over_tau_k")],
             rows=x))
    return dict(histories=rows,common_k_range_h_per_Mpc=[.015,10.],
       scope="Native adiabatic prompt-only tracer. A direct comparison of two choices of closure switch and hierarchy, not all-hierarchy convergence or an observational constraint.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("native",type=Path)
    p.add_argument("--histories",nargs="+",default=["instant","ratio5"])
    a=p.parse_args()
    print(json.dumps(run(a.native,a.histories),indent=2,allow_nan=False))
