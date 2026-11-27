"""Validate exact-mode wide RK and compare independent standard native controls."""
import argparse,hashlib,json,re
from pathlib import Path
import numpy as np
from warm_direct_vacuum_native_wide_rk_reference import MODES
from warm_light_occupation_stage_audit_reference import read


def run(folder,first):
    A=read(folder/"spec.json")
    fractions={r["label"]:r["omega_nd"]/.12 for r in A["rows"]}
    zkeys=(5,20)
    R={}
    for key,frac in fractions.items():
        l=(folder/(key+".log")).read_text()
        assert "Error" not in l and "switch on ncdm fluid" not in l
        N=float(re.search(r"Non-Cold Species Nr\..+?omega = ([^\s]+)",l)[1])
        assert abs(N/(.12*frac)-1)<6e-6
        v=[]
        for j,z in enumerate(zkeys,1):
            p=folder/f"{key}_00_z{j}_pk.dat"
            assert f"redshift z={z}" in p.read_text().splitlines()[0]
            x=np.loadtxt(p)
            i=np.searchsorted(x[:,0],MODES)
            assert np.array_equal(x[i,0],MODES)
            v.append(dict(z=z,values=x[i,1].tolist(),sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
        R[key]=v
    vector=lambda label: np.asarray([t["values"] for t in R[label]])
    out=[]
    for tier in ("compact","resolved"):
        trace=vector(tier+"-trace")
        pkey=tier+"-physical"
        akey=tier+"-amp4"
        effect=np.log(vector(pkey)/trace)
        lin=np.log(vector(akey)/trace)*(fractions[pkey]-fractions[tier+"-trace"])/(fractions[akey]-fractions[tier+"-trace"])
        out.append(dict(tier=tier,physical_minus_trace=effect.tolist(),from_amp4_tangent=lin.tolist(),
                        maximum_absolute_tangent_residual=float(np.max(abs(lin-effect)))))
    a,b=[np.asarray(t["physical_minus_trace"]) for t in out]
    compare=read(first)
    s=next(t for t in compare["trace_centered_audit"] if t["tier"]=="std")
    old=np.asarray(s["physical_minus_trace"])
    j=[compare["k_h_per_Mpc"].index(k) for k in (.3,1.,2.)]
    k=[MODES.index(k) for k in (.3,1.,2.)]
    wiggle=b[:,k]-old[:,j]
    side=np.log(vector("resolved-heldphysical")/vector("resolved-heldtrace"))-b
    assert np.all(b[:,-2:]<0.) and np.all(side[:,-2:]<0.)
    return dict(scope="Exact-mode wider early RK test of a single conditioned post-production scalar, "
                      "compares two combined hierarchy/momentum choices, not separate continuum limits.",
                native_k_h_per_Mpc=list(MODES),redshifts=list(zkeys),cases=out,
                higher_minus_lower_tier_physical_difference=(b-a).tolist(),
                std_and_wide_common_modes_h_per_Mpc=[.3,1.,2.],
                wide_resolved_minus_std_NDF_physical_difference=wiggle.tolist(),
                held_vs_direct_with_same_fraction_and_each_trace=side.tolist(),selected_raw_values=R)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("wide",type=Path);p.add_argument("std",type=Path)
    a=p.parse_args();print(json.dumps(run(a.wide,a.std),indent=2,allow_nan=False))
