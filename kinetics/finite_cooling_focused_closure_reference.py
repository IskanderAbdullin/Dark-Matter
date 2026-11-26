"""Matched precision switches for a focused native CLASS tracer audit."""
import argparse
import hashlib
import json
from pathlib import Path
import re

import numpy as np
from scipy.interpolate import PchipInterpolator

from finite_cooling_pair_split_tracer_reference import ONE, TWO, PREFIX, fields, indices, transfer

LABELS={
    "joined_on_l36":PREFIX+"_nfaon_l36_focus_sync",
    "joined_off_l36":PREFIX+"_nfaoff_focus_sync",
    "joined_off_l50":PREFIX+"_nfaoff_l50_focus_sync",
    "split_off_l36":PREFIX+"_split_nfaoff_focus_sync",
}


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def prepare(native,out):
    out.mkdir(parents=True,exist_ok=True)
    rows=[]
    for key,name in LABELS.items():
        seed=TWO if key.startswith("split") else ONE
        d=fields(native/(seed+".ini"))
        l=50 if key.endswith("50") else 36
        af=2 if "_on_" in key else 3
        d.update({"root":str(out/(name+"_")),"z_pk":"20","l_max_ncdm":str(l),
           "k_min_tau0":"24000","k_per_decade_for_pk":"5","ncdm_fluid_approximation":str(af),
           "k_output_values":"2.6944, 6.736"})
        ini=out/(name+".ini")
        ini.write_text("".join(f"{k} = {v}\n" for k,v in d.items()))
        rows.append(dict(case=key,config=ini.name,source_psd_sha256=[sha(Path(p)) for p in d["ncdm_psd_filenames"].split(",")],
          lmax=l,ncdm_fluid_approximation=af,k_wavenumber_probes_h_per_Mpc=[4.,10.]))
    return dict(cases=rows,scope="Native fixed shape and normalization. Focused accurate output uses six modes, including exact k4 and k10 h/Mpc. Lmax and absence of fluid approximation are separately recorded.")


def match(native,a,b,only):
    A,B=[fields(native/(LABELS[q]+".ini")) for q in (a,b)]
    assert set(A)==set(B)
    assert {k for k in A if A[k]!=B[k]}=={"root",only}


def compare(native,requested):
    probes=np.array([4.,10.])
    if {"joined_on_l36","joined_off_l36"}.issubset(requested):
        match(native,"joined_on_l36","joined_off_l36","ncdm_fluid_approximation")
    if {"joined_off_l36","joined_off_l50"}.issubset(requested):
        match(native,"joined_off_l36","joined_off_l50","l_max_ncdm")
    wide=native/f"{ONE}_00_z3_tk.dat"
    k,v=transfer(wide)
    v0=PchipInterpolator(np.log(k),v[:,0])(np.log(probes))
    rows=[]
    for key in requested:
        name=LABELS[key]
        log=(native/(name+".log")).read_text()
        ini=fields(native/(name+".ini"))
        assert len(re.findall("evolving mode k=",log))==6
        if "_off_" in key:
            assert "switch on ncdm fluid" not in log and ini["ncdm_fluid_approximation"]=="3"
        else:
            assert "switch on ncdm fluid" in log and ini["ncdm_fluid_approximation"]=="2"
        a,b=transfer(native/f"{name}_00_tk.dat")
        j=np.searchsorted(a,probes)
        assert np.array_equal(a[j],probes)
        x=b[j,:]
        if x.shape[1]==1:
            t=x[:,0]
            f=None
        else:
            assert key.startswith("split")
            bg=native/f"{name}_00_background.dat"
            lab=indices(bg,"rho_ncdm[0]")
            tab=np.loadtxt(bg,usecols=[lab[q] for q in ("z","(.)rho_ncdm[0]","(.)rho_ncdm[1]")])
            hot=tab[:,2]/(tab[:,1]+tab[:,2])
            f=float(np.interp(np.log(21),np.log1p(tab[::-1,0]),hot[::-1]))
            t=x@np.array([1-f,f])
        rows.append(dict(case=key,lmax=int(ini["l_max_ncdm"]),ncdm_fluid_approximation=int(ini["ncdm_fluid_approximation"]),
             native_k=a.tolist(),component_delta_over_cdm=x.tolist(),summed_tracer_delta_over_cdm=t.tolist(),
             hot_energy_fraction=f,difference_from_wide_default=(t-v0).tolist()))
    return dict(z=20,probe_k_h_per_Mpc=probes.tolist(),wide_default_reference=v0.tolist(),cases=rows,
        scope="Diagnostic native sensitivity of an extremely subdominant scalar tracer; the configuration switch can be isolated when the matching focused on/off pair is available. Sparse k and tested finite hierarchies, not total matter transfer or full physical convergence.")


def history(native):
    match(native,"joined_on_l36","joined_off_l36","ncdm_fluid_approximation")
    noms=[LABELS[k] for k in ("joined_on_l36","joined_off_l36")]
    log=(native/(noms[0]+".log")).read_text()
    out=[]
    for i,k in enumerate([2.6944,6.736]):
        filt=re.compile(r"Mode k=([.0-9e+\-]+): switch on ncdm fluid approximation at tau=([.0-9e+\-]+)")
        times=[float(tau) for kn,tau in filt.findall(log) if abs(float(kn)/k-1.)<1e-8]
        assert len(times)==1
        switch=times[0]
        controls=[]
        for n in noms:
            f=native/f"{n}_00_perturbations_k{i}_s.dat"
            with f.open() as r:
                assert abs(float(re.search(r"k = (\S+)",r.readline())[1])/k-1.)<1e-10
            names=indices(f,"delta_ncdm[0]")
            X=np.loadtxt(f,usecols=[names[q] for q in ("tau","a","delta_cdm","delta_ncdm[0]","cs2_ncdm[0]")])
            _,j=np.unique(X[:,0],return_index=True)
            X=X[j]
            assert np.all(X[:,2]<0.)
            D=X[:,3]/X[:,2]
            controls.append((X[:,0],X[:,1],D,X[:,2],D*X[:,4]))
        tau=np.unique(np.r_[np.geomspace(max(x[0][0] for x in controls)*1.001,switch*.995,240),
                       switch*np.array([.999,1.001,1.01,1.1,1.5,2.,4.,10.,30.,90.])])
        ell=np.log(tau)
        def interp(row,idx):
            x,y=np.log(row[0]),row[idx]
            keep=np.isfinite(y)
            return PchipInterpolator(x[keep],y[keep])
        D0,D1=[interp(row,2) for row in controls]
        A=interp(controls[0],1)
        C0,C1=[interp(row,3) for row in controls]
        P0,P1=[interp(row,4) for row in controls]
        diff=D1(ell)-D0(ell)
        before=float(np.max(abs(diff[tau<switch])))
        assert before<1e-5
        t20=np.exp(np.interp(np.log(1/21),np.log(controls[1][1]),np.log(controls[1][0])))
        out.append(dict(k_Mpc_minus1=k,fluid_switch_tau_Mpc=switch,scale_factor_at_switch=float(A(np.log(switch))),
            max_preswitch_ratio_difference=before,tau_probes_in_switch_units=(tau[-10:]/switch).tolist(),
            on_ratio_through_switch=D0(ell[-10:]).tolist(),off_ratio_through_switch=D1(ell[-10:]).tolist(),
            on_pressure_perturbation_over_rho_scalar_and_cdm_delta=P0(ell[-10:]).tolist(),
            off_pressure_perturbation_over_rho_scalar_and_cdm_delta=P1(ell[-10:]).tolist(),
            max_relative_cdm_change_in_switch_probes=float(np.max(np.abs(C0(ell)/C1(ell)-1.))),
            relative_cdm_change_z20=float(C0(np.log(t20))/C1(np.log(t20))-1.),
            z20_trace_on_off=[float(c(np.log(t20))) for c in (D0,D1)]))
    return dict(scope="Matched native time-resolved switch control, same scalar f0 and lmax=36; switch time read from solver logs. Not a tensor hierarchy convergence certificate.",cases=out)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("native",type=Path)
    p.add_argument("--prepare",type=Path)
    p.add_argument("--cases",choices=list(LABELS),nargs="+")
    p.add_argument("--time-probe",action="store_true")
    a=p.parse_args()
    if a.prepare is not None:
        result=prepare(a.native.resolve(),a.prepare.resolve())
    elif a.time_probe:
        result=history(a.native.resolve())
    else:
        result=compare(a.native.resolve(),a.cases)
    print(json.dumps(result,indent=2,allow_nan=False))
