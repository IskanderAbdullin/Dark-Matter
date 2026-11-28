"""Compare the exact native response under physically distinct abundance conventions."""
import argparse,json,re,hashlib
from pathlib import Path
import numpy as np
from warm_direct_vacuum_native_extended_k_reference import MODES
from warm_light_occupation_stage_audit_reference import read


def config(root,case):
    text=(root/(case+".ini")).read_text()
    r=dict((a.strip(),b.strip()) for a,b in (v.split("=",1) for v in text.splitlines() if "=" in v))
    assert r["non linear"]=="none" and r["ncdm_fluid_approximation"]=="3"
    assert r["gauge"]=="synchronous" and r["z_pk"]=="5,20"
    assert abs(float(r["omega_ncdm"])+float(r["omega_cdm"])-.12)<3e-16
    assert np.allclose(np.asarray(r["k_output_values"].split(","),float)/float(r["h"]),MODES,rtol=1e-13)
    return r


def matched(first,second):
    args={"root","omega_ncdm","omega_cdm"}
    assert {k:v for k,v in first.items() if k not in args}=={k:v for k,v in second.items() if k not in args}


def get(root,case):
    log=(root/(case+".log")).read_text()
    assert "Error" not in log and "switch on ncdm fluid" not in log and "Computing sources" in log
    c=config(root,case)
    assert "file "+c["ncdm_psd_filenames"] in log or "file native/"+Path(c["ncdm_psd_filenames"]).name in log
    omega=float(re.search(r"Non-Cold Species Nr\..+?omega = ([^\s]+)",log)[1])
    D=float(re.search(r"gives Delta N_eff = ([^\s]+)",log)[1])
    tables=[]
    for j,z in enumerate((5,20),1):
        p=root/f"{case}_00_z{j}_pk.dat"
        assert f"redshift z={z}" in p.read_text().splitlines()[0]
        x=np.loadtxt(p);ix=np.searchsorted(x[:,0],MODES)
        assert np.array_equal(x[ix,0],MODES) and np.all(x[ix,1]>0)
        tables.append(dict(z=z,native_power=x[ix,1].tolist(),sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    frac=float(c["omega_ncdm"])
    assert abs(omega/frac-1)<6e-6
    return dict(case=case,run_omega_ncdm=frac,reported_relativistic_delta_Neff=D,raw=tables)


def run(high,absolute,dominant,entropy,partial):
    records=[]
    excluded=[]
    configs={}
    for folder,category,cases in (
      (high,"held_count_rescaled",("compact-trace","compact-physical","dense-trace","dense-physical",
                            "dense-amp4","dense-heldtrace","dense-heldphysical")),
      (absolute,"absolute_vacuum_yield",("compact-absolute","dense-absolute","compact-upper_normalization")),
      (dominant,"unrelated_imposed_abundance",("ten_percent","half")),
      (entropy,"entropy_momentum_shape_only",("trace","physical","amp4"))):
        arm=[]
        for case in cases:
            if partial:
                saved=[folder/(case+".log"),folder/f"{case}_00_z1_pk.dat",folder/f"{case}_00_z2_pk.dat"]
                if not all(m.is_file() and m.stat().st_size>0 for m in saved):
                    excluded.append(dict(case=case,category=category,reason="Missing or empty native log or power table"))
                    continue
            arm.append(get(folder,case))
            configs[category,case]=config(folder,case)
        records.append(dict(category=category,cases=arm))
    raw={a["category"]:{b["case"]:b for b in a["cases"]} for a in records}
    P=lambda x:np.asarray([b["native_power"] for b in x["raw"]])
    for k in ("compact","dense"):
        cat="held_count_rescaled"
        for case in (k+"-physical",k+"-amp4"):
            if (cat,case) in configs and (cat,k+"-trace") in configs:
                matched(configs[cat,case],configs[cat,k+"-trace"])
        ab="absolute_vacuum_yield"
        for case in (k+"-absolute",k+"-upper_normalization"):
            if (ab,case) in configs and (cat,k+"-trace") in configs:
                matched(configs[ab,case],configs[cat,k+"-trace"])
    for cat,case in (("unrelated_imposed_abundance","ten_percent"),("unrelated_imposed_abundance","half")):
        if (cat,case) in configs and ("held_count_rescaled","compact-trace") in configs:
            matched(configs[cat,case],configs["held_count_rescaled","compact-trace"])
    for cat,base,cases in (("entropy_momentum_shape_only","trace",("physical","amp4")),
                             ("held_count_rescaled","dense-heldtrace",("dense-heldphysical",))):
        for case in cases:
            if (cat,case) in configs and (cat,base) in configs:
                matched(configs[cat,case],configs[cat,base])
    result={}
    for k in ("compact","dense"):
        t=raw["held_count_rescaled"]
        if {k+"-trace",k+"-physical"}<=t.keys():
            trace=P(t[k+"-trace"])
            exact=np.log(P(t[k+"-physical"])/trace)
            u=dict(rescaled_dark_composition_logP=exact.tolist())
            ab=raw["absolute_vacuum_yield"]
            if k+"-absolute" in ab:
                d1=np.log(P(ab[k+"-absolute"])/trace)
                u["absolute_direct_logP_same_CDM_fill"] = d1.tolist()
                u["ratio_of_absolute_to_rescaled_at_k10"]=(d1[:,-1]/exact[:,-1]).tolist()
            if k=="compact" and k+"-upper_normalization" in ab:
                u["upper_normalization_only_fixed_shape_logP"]=np.log(P(ab[k+"-upper_normalization"])/trace).tolist()
            result[k]=u
    check={}
    r=raw["held_count_rescaled"]
    if {"dense-heldtrace","dense-heldphysical","dense-physical","dense-trace"}<=r.keys():
        held=np.log(P(r["dense-heldphysical"])/P(r["dense-heldtrace"]))
        direct=np.array(result["dense"]["rescaled_dark_composition_logP"])
        d=r["dense-physical"]["run_omega_ncdm"]-r["dense-trace"]["run_omega_ncdm"]
        h=r["dense-heldphysical"]["run_omega_ncdm"]-r["dense-heldtrace"]["run_omega_ncdm"]
        check["dense_held_PSD_response_to_its_trace_logP"]=held.tolist()
        check["dense_held_minus_direct_slope_at_same_abundance_logP"]=(held*d/h-direct).tolist()
        check["dense_held_versus_direct_per_density_response_logP_ratio"]=(held*d/(h*direct)).tolist()
    token="absolute_direct_logP_same_CDM_fill"
    if all(token in result.get(k,{}) for k in ("compact","dense")):
        c,d=(np.array(result[k][token]) for k in ("compact","dense"))
        check["absolute_two_tiers_max_logP_difference"]=float(np.max(abs(c-d)))
        check["absolute_two_tiers_relative_response_difference"]=(c/d-1).tolist()
    if {"dense-amp4","dense-physical","dense-trace"}<=r.keys() and token in result.get("dense",{}):
        base=r["dense-trace"]["run_omega_ncdm"]
        tangent=np.log(P(r["dense-amp4"])/P(r["dense-trace"]))/(r["dense-amp4"]["run_omega_ncdm"]-base)
        actual=np.array(result["dense"][token])
        norm=raw["absolute_vacuum_yield"]["dense-absolute"]["run_omega_ncdm"]-base
        estimate=tangent*norm
        check["dense_absolute_from_tiny_abundance_tangent_logP"]=estimate.tolist()
        check["dense_absolute_over_tiny_abundance_tangent"]=(actual/estimate).tolist()
    dom={}
    t=raw["held_count_rescaled"]
    if "compact-trace" in t:
        for k,a in raw["unrelated_imposed_abundance"].items():
            p=P(a)/P(t["compact-trace"])
            dom[k]=dict(power_ratios=p.tolist(),sampled_half_level_brackets=[
                         [[MODES[i],MODES[i+1]] for i in range(len(MODES)-1) if (v[i]-.5)*(v[i+1]-.5)<0]
                         for v in p])
    e=raw["entropy_momentum_shape_only"]
    cooling={}
    if {"trace","physical"}<=e.keys():
        c=np.log(P(e["physical"])/P(e["trace"]))
        cooling["cooling_physical_minus_own_trace"]=c.tolist()
        if "compact" in result:
            t=np.asarray(result["compact"]["rescaled_dark_composition_logP"])
            cooling["cooler_over_original_at_k10"]= (c[:,-1]/t[:,-1]).tolist()
        if "amp4" in e:
            a=np.log(P(e["amp4"])/P(e["trace"]))
            frac=(e["physical"]["run_omega_ncdm"]-e["trace"]["run_omega_ncdm"])/(e["amp4"]["run_omega_ncdm"]-e["trace"]["run_omega_ncdm"])
            cooling["max_tangent_difference"]=float(max(abs(a*frac-c).ravel()))
    return dict(scope="Exact modes in linear native CLASS, distinct abundance and momentum-shape hypotheses. "
                      "The absolute direct case independently fills missing stable abundance with CDM; "
                      "no such component is generated or fitted by this calculation.",
                k_h_per_Mpc=list(MODES),z=[5.,20.],partial=partial,
                raw_cases=records,excluded_incomplete=excluded,matched_comparisons=result,numerical_cross_checks=check,imposed_counterfactual=dom,entropy_counterfactual=cooling)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("high","absolute","dominant","entropy"):
        p.add_argument(k,type=Path)
    p.add_argument("--partial",action="store_true")
    a=p.parse_args();print(json.dumps(run(a.high,a.absolute,a.dominant,a.entropy,a.partial),indent=2,allow_nan=False))
