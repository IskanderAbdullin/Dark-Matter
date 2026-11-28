"""Keep the fixed-metric contrast separate from the spectrum smoothing error."""
import argparse
import json
from pathlib import Path

import numpy as np

from finite_cooling_matched_history_los_reference import KEY, sha


def run(eos, pair_export, prompt_export, narrow_recoil, wide_recoil, wide_paths, narrow_paths):
    read=lambda p: json.loads(p.read_text())
    mix,prompt,nr,wr=map(read,(pair_export,prompt_export,narrow_recoil,wide_recoil))
    assert nr["sigma_ln_q"]==.10 and nr["source_sha256"]==wr["source_sha256"]==mix["pairs_sha256"]
    assert sha(eos)==nr["eos_sha256"]==wr["eos_sha256"]==mix["eos_sha256"]==prompt["eos_sha256"]
    rows={}
    for c,a,b in zip(("instant","ratio5"),wide_paths,narrow_paths,strict=True):
        W,N=read(a),read(b)
        assert W["z"]==N["z"]==20.
        assert nr["components"][c]["rest_omega_h2"]==wr["components"][c]["rest_omega_h2"]
        for D,R,sigma in ((W,wr,.14),(N,nr,.10)):
            j=[x for x in mix["components"] if x["case"]==c and x["sigma_ln_q"]==sigma]
            pr=[x for x in prompt["components"] if x["case"]==c and x["sigma_ln_q"]==sigma and x["species"]=="scalar"]
            assert len(j)==len(pr)==1
            assert D["spectrum_sha256"]==[j[0]["sha256"],pr[0]["sha256"],R["components"][c]["table_sha256"]]
        xW={round(r["k_h_per_Mpc"],12):r for r in W["results"]}
        assert list(np.round([r["k_h_per_Mpc"] for r in N["results"]],12))==[4.,10.,22.]
        cases=[]
        for Y in N["results"]:
            k=round(Y["k_h_per_Mpc"],12)
            X=xW[k]
            assert all(X[key]==Y[key] for key in ("metric_sha256","original_native_metric_sha256","tau_final_Mpc","initial_k_tau"))
            w,n=[float(D["cases"][-1][KEY][0][0]) for D in (X,Y)]
            cases.append(dict(k_h_per_Mpc=k,joined_wide=w,joined_narrow=n,change_narrow_minus_wide=n-w,
                fractional_change_narrow_minus_wide=n/w-1.,shared_metric_sha256=X["original_native_metric_sha256"],
                max_numerical_refinement=max(D["max_time_and_momentum_refinement_difference"] for D in (X,Y))))
        rows[c]=cases
    a,b=[[z["joined_narrow"] for z in rows[c]] for c in ("instant","ratio5")]
    delta=(np.array(b)-a).tolist()
    assert min(delta)>.5
    return dict(rows=rows,narrow_finite_minus_instant=delta,raw_record_sha256=[sha(p) for p in [*wide_paths,*narrow_paths]],
                scope="Width stress test for positive moment-controlled smoothing (0.14 vs 0.10) at three identical fixed-metric probes. Not continuum source accuracy or a total-matter effect.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for name in ("eos","pair_export","prompt_export","narrow_recoil","wide_recoil"):
        p.add_argument(name,type=Path)
    p.add_argument("--wide",nargs=2,type=Path,required=True)
    p.add_argument("--narrow",nargs=2,type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.pair_export,a.prompt_export,a.narrow_recoil,a.wide_recoil,a.wide,a.narrow),indent=2,allow_nan=False))
