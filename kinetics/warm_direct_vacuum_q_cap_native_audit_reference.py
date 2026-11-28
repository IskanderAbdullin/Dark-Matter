"""Audit same-step tail shortening and higher radial q sampling on exact native modes."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from warm_light_occupation_stage_audit_reference import read


def pair(folder,names,K):
    vals=[]
    for j,z in enumerate((5,20),1):
        rows=[];sha=[]
        for l in names:
            f=folder/f"{l}_00_z{j}_pk.dat"
            assert f"redshift z={z}" in f.read_text().splitlines()[0]
            a=np.loadtxt(f)
            I=np.searchsorted(a[:,0],K)
            assert np.array_equal(a[I,0],K)
            assert np.all(a[I,1]>0)
            rows.append(a[I,1])
            sha.append(hashlib.sha256(f.read_bytes()).hexdigest())
        vals.append(dict(z=z,physical_minus_trace=np.log(rows[1]/rows[0]).tolist(),
                         exact_mode_trace_and_physical_sha256=sha))
    return vals


def run(wide,fine,cap):
    W=read(wide)
    K=W["native_k_h_per_Mpc"]
    controls=[]
    Q=read(cap/"spec.json")["cases"]
    new=pair(fine,("trace","physical"),K)
    ns=np.asarray([p["physical_minus_trace"] for p in new])
    old=np.asarray(next(a for a in W["cases"] if a["tier"]=="compact")["physical_minus_trace"])
    q96=np.asarray(next(a for a in W["cases"] if a["tier"]=="resolved")["physical_minus_trace"])
    for r,full in (("step70",old),("step168",ns)):
        cut=pair(cap,(r+"-trace",r+"-physical"),K)
        short=np.asarray([a["physical_minus_trace"] for a in cut])
        q=next(a for a in Q if a["tag"]==r+"-trace")
        passed=bool(np.max(abs(short-full))<1e-8)
        controls.append(dict(sampling=r,initial_all_mode_1e_minus8_target_met=passed,retained_q_max=q["cut_q"],unchanged_linear_step=q["mesh_step"],
            retained_number_and_q2_fraction=q["retained_number_and_q2_fraction"],
            max_exact_power_response_change_due_only_to_shorter_q_cap=float(np.max(abs(short-full))),
            changes_due_to_shorter_q_cap=(short-full).tolist(),shorter_q_results=cut))
    return dict(scope="Early exact-mode RK values, separating closer radial spacing from removal of a negligibly "
                      "populated high-q region. Native models renormalize to the same present scalar abundance.",
                k_h_per_Mpc=K,full_cap_168_results=new,cap_controls=controls,
                full_168_minus_96=(ns-q96).tolist())


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("wide","full168","cap"):
        p.add_argument(k,type=Path)
    a=p.parse_args();print(json.dumps(run(a.wide,a.full168,a.cap),indent=2,allow_nan=False))
