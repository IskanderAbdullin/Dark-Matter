"""Native conditional CLASS response to reweighting an archived scalar shape."""
import argparse
import json
from pathlib import Path
import re

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

from archived_class_pressure_balanced_reference import read,at

FRACTIONS={"000":0.,"0002":.002,"001":.01,"010":.10,"050":.50,"all":1.}
OMEGA_B=.02237


def case_name(a,raw=False):
    return f"scalar_fraction_{a}"+("_raw" if raw else "_smooth")


def prepare(out,fine):
    from archived_class_scalar_tracer_reference import generate
    values=generate(fine,out)
    om=values["dark_density_omega_h2"]
    base=(out/"scalar_tracer_smooth128_wide_sync.ini").read_text()
    ref="scalar_tracer_smooth128_wide_sync"
    paths=[]
    for a,f in FRACTIONS.items():
        for raw in (False,True) if a in ("001","010") else (False,):
            name=case_name(a,raw)
            s=base.replace(ref+"_",name+"_")
            c=max(om*(1-f),1e-8)
            if f==0.:
                s=s[:s.index("omega_cdm =")]+f"omega_cdm = {om:.17g}\nk_min_tau0 = 100\n"
            else:
                s=re.sub(r"omega_cdm = [^\n]+",f"omega_cdm = {c:.17g}",s)
                s=re.sub(r"omega_ncdm = [^\n]+",f"omega_ncdm = {om-c:.17g}",s)
            if raw:
                s=re.sub(r"ncdm_psd_filenames = [^\n]+",f"ncdm_psd_filenames = {fine/'scalar.dat'}",s)
                s=s.replace("ncdm_N_momentum_bins = 128","ncdm_N_momentum_bins = 160")
                s=s.replace("ncdm_maximum_q = 96","ncdm_maximum_q = 63.8376753507014")
                s=s.replace("l_max_ncdm = 10","l_max_ncdm = 14")
                s=s.replace("ncdm_fluid_trigger_tau_over_tau_k = 12","ncdm_fluid_trigger_tau_over_tau_k = 20")
            p=out/(name+".ini")
            p.write_text(s)
            paths.append(str(p))
    return dict(generated_inputs=paths,dark_density_omega_h2=om)


def transfer(out,case,i,z,k):
    p=out/f"{case}_00_z{i}_tk.dat"
    with p.open() as f:
        assert re.search(rf"redshift z={z}(?:[\s.]|$)",f.readline())
        labels=next(x for x in f if x.startswith("#") and "d_b" in x and "d_cdm" in x)
    names={j:int(i)-1 for i,j in re.findall(r'(\d+):(\S+)',labels)}
    a=np.loadtxt(p)
    assert a[0,0]<=min(k)<=max(k)<=a[-1,0]
    v=PchipInterpolator(np.log(a[:,0]),a[:,1:],axis=0)(np.log(k))
    return [v[:,names[key]-1] if key in names else None for key in ("d_b","d_cdm","d_ncdm[0]")]


def analyze(out,codes,raw_codes,fine):
    p=json.loads((fine/"manifest.json").read_text())
    om=p["components"][0]["omega_h2_from_retained_rest_count"]+1e-8
    pr=np.array([.1,.3,1.,2.,4.,8.,10.])
    k=np.unique(np.r_[np.geomspace(.03,12.,90),pr])
    q=np.searchsorted(k,pr)
    p0=case_name("000")
    rows=[]
    choices=[(a,False) for a in codes]+[(a,True) for a in raw_codes]
    for i,z in enumerate((0,5,20),1):
        A=at(read(out,p0,i),k)
        b0,c0,_=transfer(out,p0,i,z,k)
        total0=OMEGA_B*b0+om*c0
        assert np.all(total0<0)
        d={}
        for a,israw in choices:
            n=case_name(a,israw)
            f=FRACTIONS[a]
            cd=max(om*(1.-f),1e-8)
            frac=(om-cd)/om
            P=at(read(out,n,i),k)/A
            b,c,s=transfer(out,n,i,z,k)
            cold=OMEGA_B*b+om*c
            growth=cold/total0
            weight=1.+(om-cd)*(s-c)/cold
            signed=growth*weight
            half=np.where((P[:-1]>.5)&(P[1:]<=.5))[0]
            p_half=None if not len(half) else float(brentq(
                lambda x: float(at(read(out,n,i),np.array([x]))[0]/at(read(out,p0,i),np.array([x]))[0])-.5,
                k[half[0]],k[half[0]+1]))
            d[n]=dict(actual_scalar_fraction=frac,
                total_linear_power_ratio=P[q].tolist(),
                change_in_cold_and_baryon_growth_amplitude=growth[q].tolist(),
                direct_matter_weighting_amplitude=weight[q].tolist(),
                synchronous_weighted_vs_native_power_ratio_max_error_probes_above_0p05=
                    float(max(abs(((signed[q]**2)/P[q]-1)[P[q]>.05]),default=0.)),
                probes_below_0p05_count=int(np.count_nonzero(P[q]<=.05)),
                first_bracketed_power_half_h_per_Mpc=p_half,
                dense_power_ratio=[float("%.8g"%v) for v in P])
        small=[case_name(a) for a in ("0002","001") if case_name(a) in d]
        slopes={n:(np.log(np.asarray(d[n]["total_linear_power_ratio"]))/d[n]["actual_scalar_fraction"]).tolist()
                for n in small}
        rows.append(dict(z=z,k_h_per_Mpc=pr.tolist(),dense_k_grid_h_per_Mpc=[float("%.8g"%v) for v in k],
                         finite_small_fraction_log_power_slopes=slopes,cases=d))
    return dict(omega_dark=om,omega_b=OMEGA_B,rows=rows,
       scope="Counterfactual frozen scalar shape, fixed total dark density and adiabatic CLASS histories. Each fraction changes the gravitational background and growth. None of these fractions is the actual low-reheat production result (fraction about 2.53e-9). Growth/weight decomposition uses the synchronous late matter contrasts, not a fit or a Ly-alpha comparison.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("out",type=Path)
    p.add_argument("fine",type=Path)
    p.add_argument("--prepare",action="store_true")
    p.add_argument("--codes",nargs="+",default=["001","010","050","all"])
    p.add_argument("--raw-codes",nargs="*",default=[])
    a=p.parse_args()
    out=a.out.resolve()
    out.mkdir(exist_ok=True,parents=True)
    result=prepare(out,a.fine.resolve()) if a.prepare else analyze(out,a.codes,a.raw_codes,a.fine)
    print(json.dumps(result,indent=2,allow_nan=False))
