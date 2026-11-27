"""Resolve individual momentum, multipole and interaction shifts of the wide native tests."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from warm_light_occupation_stage_audit_reference import read


def run(wide,cross):
    root=read(wide)
    K=root["native_k_h_per_Mpc"]
    E={r["tier"]:np.asarray(r["physical_minus_trace"]) for r in root["cases"]}
    extra={}
    for k in ("q_only","l_only"):
        res=[]
        for j,z in enumerate(root["redshifts"],1):
            allp=[]
            sources=[]
            for arm in ("trace","physical"):
                tag=k+"-"+arm
                name=cross/f"{tag}_00_z{j}_pk.dat"
                a=np.loadtxt(name)
                i=np.searchsorted(a[:,0],K)
                assert np.array_equal(a[i,0],K) and np.all(a[i,1]>0.)
                assert f"redshift z={z}" in name.read_text().splitlines()[0]
                allp.append(a[i,1]);sources.append(hashlib.sha256(name.read_bytes()).hexdigest())
                assert "Error" not in (cross/(tag+".log")).read_text()
            res.append(dict(z=z,physical_minus_trace=np.log(allp[1]/allp[0]).tolist(),
                            trace_and_physical_sha256=sources))
        extra[k]=res
    A,B=E["compact"],E["resolved"]
    cq=np.asarray([r["physical_minus_trace"] for r in extra["q_only"]])
    cl=np.asarray([r["physical_minus_trace"] for r in extra["l_only"]])
    mix=B-cq-cl+A
    assert np.all(np.stack((A,B,cq,cl))[:,:,-2:]<0.)
    return dict(scope="A finite two-axis resolution cross at identical exact modes, same early start, no fluid closure "
                      "or gauge change. The interaction term is not a continuum error estimate.",
                k_h_per_Mpc=K,redshifts=root["redshifts"],other_corners=extra,
                shifts_from_compact=dict(q_only=(cq-A).tolist(),l_only=(cl-A).tolist(),
                                         resolved=(B-A).tolist(),non_additive=mix.tolist()))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("wide",type=Path);p.add_argument("cross",type=Path)
    a=p.parse_args();print(json.dumps(run(a.wide,a.cross),indent=2,allow_nan=False))
