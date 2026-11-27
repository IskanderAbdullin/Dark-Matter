"""Matched-density versus own absolute-density exact modes for a reheating ladder."""
import argparse
import json
from pathlib import Path

import numpy as np

from warm_direct_vacuum_native_extended_k_reference import MODES
from warm_direct_vacuum_native_normalization_audit_reference import get,config,matched
from warm_light_occupation_stage_audit_reference import read


def run(spec,root,requested):
    records=[]
    response={}
    b=read(spec)
    powers=lambda r:np.asarray([j["native_power"] for j in r["raw"]])
    expected50=sum(next(x for x in b["rows"] if x["T_reheat_GeV"]==.05)["omega_h2_per_source"])
    for row in b["rows"]:
        t=row["T_reheat_GeV"]*1000
        tag=f"direct{t:.0f}"
        if int(round(t)) not in requested:
            continue
        names=("trace","absolute","matched50") if t!=50. else ("trace","absolute")
        vals={n:get(root,tag+"_"+n) for n in names}
        base=config(root,tag+"_trace")
        for name in names:
            matched(base,config(root,tag+"_"+name))
        absolute=powers(vals["absolute"])
        trace=powers(vals["trace"])
        d0=np.log(absolute/trace)
        amp=vals["absolute"]["run_omega_ncdm"]
        np.testing.assert_allclose(amp,sum(row["omega_h2_per_source"]),rtol=6e-14,atol=0)
        d1=np.log(powers(vals["matched50"])/trace) if "matched50" in vals else d0.copy()
        if "matched50" in vals:
            assert abs(vals["matched50"]["run_omega_ncdm"]/expected50-1.)<3e-6
        response[tag]=(d0,d1,trace)
        records.append(dict(T_reheat_MeV=t,raw_runs=list(vals.values()),
                            density_relative_offset_from_native_50_for_matched_case=(vals["matched50"]["run_omega_ncdm"]/expected50-1) if "matched50" in vals else 0.,
                            log_power_change_at_own_abundance=d0.tolist(),
                            log_power_change_at_matched_50_MeV_abundance=d1.tolist()))
    additions=[]
    if "direct50" in response:
        _,base,tr=response["direct50"]
        for tag,(a,b,trace) in response.items():
            shape=b-base
            norm=a-b
            both=a-base
            np.testing.assert_allclose(shape+norm,both,atol=1e-16)
            additions.append(dict(tag=tag,shape_only_logP_at_matched_50_abundance_vs_T50=shape.tolist(),
                        additional_normalization_logP_at_own_shape=norm.tolist(),
                        sum_change_in_logP_after_recentring_on_own_traces=both.tolist(),
                        trace_power_difference_logP_vs_T50=np.log(trace/tr).tolist()))
    assert len(records)==len(requested)
    return dict(k_h_per_Mpc=MODES,z=[5,20],cases=records,decomposition=additions,
                scope="Conditional linear mPk with total omega_dark=.12. Equal-abundance subtraction "
                      "compares the shapes at the 50 MeV scalar density, up to the recorded central quadrature offset. Subsequent change "
                      "uses each shape's own direct scalar yield and each trace includes the rest as CDM.")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("spec",type=Path)
    p.add_argument("root",type=Path)
    p.add_argument("temps",nargs="+",type=int)
    a=p.parse_args()
    print(json.dumps(run(a.spec,a.root,set(a.temps)),indent=2,allow_nan=False))
