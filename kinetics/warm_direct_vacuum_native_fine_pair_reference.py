"""Audit the completed full fine NDF physical/trace pair separately from queued amplifiers."""
import argparse,hashlib,json,re
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator,CubicSpline
from warm_light_occupation_stage_audit_reference import read


def run(folder,wide,prior):
    spec=read(folder/"spec.json")
    E={a["label"]:a for a in spec["cases"]}
    for k in ("fine-physical","fine-trace"):
        t=E[k]
        assert t["quadrature"]==128 and t["lmax"]==32 and t["rtol"]==8e-10
        log=(folder/(k+".log")).read_text()
        assert "Error" not in log and "switch on ncdm fluid" not in log
        rho=float(re.search(r"Non-Cold Species Nr\..+?omega = ([^\s]+)",log)[1])
        assert abs(rho/(.12*t["rest_fraction"])-1)<6e-6
    W=read(wide)
    H=read(prior)
    K=np.asarray(H["k_h_per_Mpc"])
    x=np.log(K)
    upper=np.asarray(next(r for r in W["cases"] if r["tier"]=="resolved")["physical_minus_trace"])
    old=np.asarray(next(r for r in H["trace_centered_audit"] if r["tier"]=="std")["physical_minus_trace"])
    ix=np.searchsorted(K,[.3,1.,2.])
    jx=[W["native_k_h_per_Mpc"].index(k) for k in (.3,1.,2.)]
    rows=[]
    for j,z in enumerate((5,20),1):
        R=[]
        for key in ("fine-trace","fine-physical"):
            p=folder/f"{key}_00_z{j}_pk.dat"
            assert f"redshift z={z}" in p.read_text().splitlines()[0]
            a=np.log(np.loadtxt(p))
            assert np.all(np.diff(a[:,0])>0.)
            R.append((PchipInterpolator(a[:,0],a[:,1])(x),CubicSpline(a[:,0],a[:,1])(x),
                      hashlib.sha256(p.read_bytes()).hexdigest()))
        signal=R[1][0]-R[0][0]
        alternate=R[1][1]-R[0][1]
        rows.append(dict(z=z,fine_trace_physical_sha256=[a[2] for a in R],
             physical_minus_trace=signal.tolist(),cubic_minus_pchip=(alternate-signal).tolist(),
             fine_minus_wide_RK_common_0p3_1_2=(signal[ix]-upper[j-1,jx]).tolist(),
             fine_minus_std_NDF=(signal-old[j-1]).tolist()))
    return dict(scope="Full fine NDF trace and physical pair only; not a claim that the higher-order amplified "
                      "or held-shape histories are complete. Fixed conditional model.",
                native_requested_comparison_k_h_per_Mpc=K.tolist(),cases=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("native","wide","prior"):
        p.add_argument(k,type=Path)
    a=p.parse_args();print(json.dumps(run(a.native,a.wide,a.prior),indent=2,allow_nan=False))
