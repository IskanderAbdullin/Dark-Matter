"""Native scalar-tracer configurations at the fixed-count benchmark abundances."""
import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np

from archived_class_scalar_tracer_reference import extract
from archived_collisionless_streaming_reference import Background
from class_exact_coordinate_remap_reference import scan


def prepare(psd,manifest,out):
    out.mkdir(parents=True,exist_ok=True)
    source=json.loads(manifest.read_text())
    rows=[]
    for a in source["components"]:
        if a["species"]!="scalar":
            continue
        key,sigma=a["case"],a["sigma_ln_q"]
        p=psd/a["file"]
        assert hashlib.sha256(p.read_bytes()).hexdigest()==a["sha256"]
        n=160 if sigma==.10 else 224
        ranking=scan(p,n)
        g=ranking["at_native_scale_1"] if key=="instant" else ranking["best"]
        assert g["maximum_residual"]<1e-6
        scale=g["scale"]
        name=f"finite_cooling_tracer_{key}_s{round(sigma*100)}_sync"
        dst=out/(name+".dat")
        dat=np.loadtxt(p)
        dat[:,0]/=scale
        np.savetxt(dst,dat,fmt="%.16e")
        abundance=a["actual_fixed_count_Omega_h2"]
        assert 0<abundance<2.4e-9
        (out/(name+".ini")).write_text(f"""output = mPk,dTk
P_k_max_h/Mpc = 10
z_pk = 0,5,20
h = 0.6736
T_cmb = 2.72548
omega_b = 0.02237
A_s = 2.1e-9
n_s = 0.965
tau_reio = 0.054
N_ur = 3.044
gauge = synchronous
background_verbose = 2
thermodynamics_verbose = 0
perturbations_verbose = 3
write background = yes
headers = yes
k_per_decade_for_pk = 14
tol_perturbations_integration = 1e-07
root = {out}/{name}_
omega_cdm = 0.12
N_ncdm = 1
m_ncdm = 86
T_ncdm = {scale}
deg_ncdm = 1
use_ncdm_psd_files = 1
ncdm_psd_filenames = {dst}
tol_ncdm = 1e-05
tol_ncdm_bg = 1e-07
omega_ncdm = {abundance}
ncdm_quadrature_strategy = 2
ncdm_maximum_q = 20
ncdm_N_momentum_bins = {n}
ncdm_fluid_trigger_tau_over_tau_k = 20
l_max_ncdm = 14
k_min_tau0 = 100
""")
        rows.append(dict(name=name,case=key,sigma=sigma,native_momentum_nodes=n,
            input_sha256=a["sha256"],remap_sha256=hashlib.sha256(dst.read_bytes()).hexdigest(),
            momentum_scale=scale,fixed_omega_h2=abundance,external_moment_gate=g))
    return dict(cases=rows,scope="Adiabatic collisionless scalar tracer after 1 MeV. The dominant 20 MeV species is replaced by a cold component, not evolved jointly; reconstructed smoothing and external quadrature checks do not certify perturbations.")


def compare(out,eos,cases):
    ai=float(Background(eos).scale(.001))
    v=extract(out,cases,ai)
    for r in v["anchor_diagnostics"]:
        line=(out/(r["case"]+".ini")).read_text()
        a=float(re.search(r"(?m)^omega_ncdm = ([^\n]+)$",line)[1])
        e=abs(r["present_NATIVE_energy_density_omega_h2"]/a-1.)
        assert e<2e-9
        assert r["largest_mode_k_times_tau_at_archive"]<.001
        r["target_fixed_count_scalar_omega_h2"]=a
        r["relative_native_density_defect"]=e
    for row in v["rows"]:
        rows=[]
        for tag in ("instant","ratio5"):
            picks=[c for c in row["cases"] if ("_"+tag+"_") in c["case"]]
            if not picks:
                continue
            lo,hi=max(c["native_k_h_per_Mpc"][0] for c in picks),min(c["native_k_h_per_Mpc"][-1] for c in picks)
            assert lo<.015 and hi>10.
            from scipy.interpolate import PchipInterpolator
            grid=np.geomspace(.015,10.,201)
            q=np.array([PchipInterpolator(np.log(c["native_k_h_per_Mpc"]),
                                c["native_delta_scalar_over_delta_cdm"])(np.log(grid)) for c in picks])
            e=float(np.max(np.ptp(q,axis=0)))
            if len(picks)==2:
                assert e<.016
            if tag=="ratio5":
                assert np.min(q)>.5
            rows.append(dict(history=tag,smoothing_control_count=len(picks),max_spread_between_shape_and_momentum_grid_controls=e,
                    common_test_range_h_per_Mpc=[.015,10.],minimum_amplitude_on_control_grid=float(q.min())))
        row["per_history_control_comparison"]=rows
        del row["probe_max_variant_spread"]
        for c in row["cases"]:
            c["native_k_max_h_per_Mpc"]=c["native_k_h_per_Mpc"][-1]
            del c["native_k_h_per_Mpc"],c["native_delta_scalar_over_delta_cdm"]
    v["scope"]="Adiabatic collisionless scalar tracer: scalar contrast relative to CDM in the same run, cold approximation for dominant stable matter, actual prompt-only scalar density. Not total-matter suppression, pre-1MeV perturbations, a dominant scalar transfer, or an observational constraint."
    return v


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("out",type=Path)
    p.add_argument("--psd",type=Path)
    p.add_argument("--manifest",type=Path)
    p.add_argument("--eos",type=Path)
    p.add_argument("--compare",nargs="+")
    a=p.parse_args()
    if a.compare:
        assert a.eos is not None
        result=compare(a.out.resolve(),a.eos,a.compare)
    else:
        assert a.psd is not None and a.manifest is not None
        result=prepare(a.psd.resolve(),a.manifest,a.out.resolve())
    print(json.dumps(result,indent=2,allow_nan=False))
