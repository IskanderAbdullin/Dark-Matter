"""Density-weighted fixed-metric test-particle contrast, not a transfer solve."""
import argparse
import json
from pathlib import Path

import numpy as np

from finite_cooling_matched_history_los_reference import KEY, sha


def read(p):
    return json.loads(p.read_text())


def run(stable,double,input_catalog,spectrum,retarget,driver,scalar):
    C,J,A,F,R,V=map(read,(stable,double,input_catalog,spectrum,retarget,driver))
    assert C["eos_sha256"]==J["eos_sha256"]==A["eos_sha256"]==F["eos_sha256"]
    histories=("instant","ratio5")
    assert C["spectrum_sha256"]==J["spectrum_sha256"]==[
        next(z["sha256"] for z in A["components"] if (z["case"],z["species"],z["sigma_ln_q"])==(x,"stable",.14))
        for x in histories]
    opts=V["native_parameters"]
    dark,baryon=float(opts["omega_cdm"]),float(opts["omega_b"])
    assert abs(dark-.120)<1e-12
    out=[]
    for i,(case,raw,rc) in enumerate(zip(histories,scalar,R["cases"],strict=True)):
        A=read(raw)
        assert A["z"]==20. and A["eos_sha256"]==C["eos_sha256"]
        M=next(x for x in F["components"] if x["case"]==case and x["sigma_ln_q"]==.14)
        assert A["spectrum_sha256"][0]==M["sha256"]
        fs=M["actual_fixed_count_Omega_h2"]/(rc["recovered_stable_omega"]+M["actual_fixed_count_Omega_h2"])
        rows=[]
        for n,t in zip(C["results"],A["results"],strict=True):
            k=n["k_h_Mpc"]
            assert abs(k-t["k_h_per_Mpc"])<1e-12
            assert n["metric_sha256"]==t["original_native_metric_sha256"]
            s=t["cases"][-1][KEY][0][0]
            d=[]
            for v in n["data"][-2:]:
                assert v["mass_eV"]==2e7
                dark_amp=v["cold_delta"]
                base=dark*dark_amp+baryon*v["exact_baryon_delta"]
                c=(dark*dark_amp)/base
                ratio1=(1.-fs)*v["kinetic_minus_cold_over_cold"][i]
                ratio2=fs*(s/v["cold_over_native"]-1.)
                rel=c*(ratio1+ratio2)
                # Check arithmetic against independent, less stable direct summation.
                pert=dark*dark_amp*((1.-fs)*(1.+v["kinetic_minus_cold_over_cold"][i])+fs*s/v["cold_over_native"])
                alt=((pert+baryon*v["exact_baryon_delta"])/base)**2-1.
                assert abs(alt-rel*(2.+rel))<2e-14
                d.append(dict(start_k_tau=v["start_k_tau"],stable_density_contrast_correction=ratio1,
                              scalar_density_contrast_correction=ratio2,baryon_dilution=c,
                              direct_fixed_metric_power_ratio_minus_one=rel*(2.+rel)))
            rows.append(dict(k_h_Mpc=k,probes=d))
        out.append(dict(history=case,scalar_fraction_of_fixed_dark_density=fs,rows=rows))
    d=J["results"][-1]
    g=C["results"][-1]
    assert d["metric_sha256"]==g["metric_sha256"] and d["k_h_Mpc"]==g["k_h_Mpc"]==22.
    X,Y=[np.array([r["kinetic_minus_cold_over_cold"] for r in z["data"][-2:]]) for z in (g,d)]
    mass_defect=float(np.max(abs((X-4*Y))))
    assert mass_defect<2e-12
    return dict(data=out,control_20_vs_40_MeV_max_abs_m_inverse2_defect=mass_defect,
         sha256={p.name:sha(p) for p in (stable,double,input_catalog,spectrum,retarget,driver,*scalar)},
         scope="Direct linear density contrast composition with the same external metric, native baryon perturbation and dark composition weights. The cold comparator is numerical on that metric; no Einstein backreaction, self-consistent matter transfer, CMB, or data constraints.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for k in ("stable","double","input_catalog","spectrum","retarget","driver"):
        p.add_argument(k,type=Path)
    p.add_argument("scalar",type=Path,nargs=2)
    a=p.parse_args()
    print(json.dumps(run(a.stable,a.double,a.input_catalog,a.spectrum,a.retarget,a.driver,a.scalar),indent=2,allow_nan=False))
