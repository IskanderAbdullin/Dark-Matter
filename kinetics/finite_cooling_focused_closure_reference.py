"""Matched precision switches for a focused native CLASS tracer audit."""
import argparse
import hashlib
import json
from pathlib import Path
import re

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

from finite_cooling_pair_split_tracer_reference import ONE, TWO, PREFIX, fields, indices, transfer

LABELS={
    "joined_on_l36":PREFIX+"_nfaon_l36_focus_sync",
    "joined_off_l36":PREFIX+"_nfaoff_focus_sync",
    "joined_off_l50":PREFIX+"_nfaoff_l50_focus_sync",
    "split_off_l36":PREFIX+"_split_nfaoff_focus_sync",
    "split_off_l36_rk":PREFIX+"_split_nfaoff_l36_rk_focus_sync",
    "duplicate_on_l36":PREFIX+"_duplicate_on_l36_focus_sync",
    "duplicate_on_l36_rk":PREFIX+"_duplicate_on_l36_rk_focus_sync",
    "joined_on_l36_rk":PREFIX+"_nfaon_l36_rk_focus_sync",
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
        if key.startswith("duplicate"):
            d["N_ncdm"]="2"
            w=float(d["omega_ncdm"])/2
            d["omega_ncdm"]=f"{w:.17g},{w:.17g}"
            for k in ("m_ncdm","T_ncdm","deg_ncdm","use_ncdm_psd_files","ncdm_psd_filenames",
                      "ncdm_quadrature_strategy","ncdm_maximum_q","ncdm_N_momentum_bins"):
                d[k]=d[k]+","+d[k]
        if key.endswith("_rk"):
            d["evolver"]="0"
        ini=out/(name+".ini")
        ini.write_text("".join(f"{k} = {v}\n" for k,v in d.items()))
        rows.append(dict(case=key,config=ini.name,source_psd_sha256=[sha(Path(p)) for p in d["ncdm_psd_filenames"].split(",")],
          lmax=l,ncdm_fluid_approximation=af,explicit_rk=key.endswith("_rk"),k_wavenumber_probes_h_per_Mpc=[4.,10.]))
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
    if {"joined_on_l36_rk","duplicate_on_l36_rk"}.issubset(requested):
        A,B=[fields(native/(LABELS[q]+".ini")) for q in ("joined_on_l36_rk","duplicate_on_l36_rk")]
        names=("m_ncdm","T_ncdm","deg_ncdm","use_ncdm_psd_files","ncdm_psd_filenames",
               "ncdm_quadrature_strategy","ncdm_maximum_q","ncdm_N_momentum_bins")
        assert A.keys()==B.keys() and A["N_ncdm"]=="1" and B["N_ncdm"]=="2"
        assert all(A[k]==B[k] for k in A if k not in (*names,"root","N_ncdm","omega_ncdm"))
        assert all([A[k],A[k]]==B[k].split(",") for k in names)
        assert all(abs(2*float(y)/float(A["omega_ncdm"])-1.)<1e-14 for y in B["omega_ncdm"].split(","))
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
            assert key.startswith(("split","duplicate"))
            bg=native/f"{name}_00_background.dat"
            lab=indices(bg,"rho_ncdm[0]")
            tab=np.loadtxt(bg,usecols=[lab[q] for q in ("z","(.)rho_ncdm[0]","(.)rho_ncdm[1]")])
            hot=tab[:,2]/(tab[:,1]+tab[:,2])
            f=float(np.interp(np.log(21),np.log1p(tab[::-1,0]),hot[::-1]))
            t=x@np.array([1-f,f])
        rows.append(dict(case=key,lmax=int(ini["l_max_ncdm"]),ncdm_fluid_approximation=int(ini["ncdm_fluid_approximation"]),
             native_k=a.tolist(),component_delta_over_cdm=x.tolist(),summed_tracer_delta_over_cdm=t.tolist(),
             hot_energy_fraction=f,explicit_rk=ini.get("evolver")=="0",difference_from_wide_default=(t-v0).tolist()))
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


def wider_default(native):
    B=PREFIX+"_wide_higher_sync"
    a,b=[fields(native/(n+".ini")) for n in (ONE,B)]
    D={k:[a[k],b[k]] for k in a if a[k]!=b[k] and k!="root"}
    assert set(D)=={"P_k_max_h/Mpc","l_max_ncdm","ncdm_fluid_trigger_tau_over_tau_k"}
    rows=[]
    for i,z in enumerate((0,5,20),1):
        name=lambda c:native/f"{c}_00_z{i}_tk.dat"
        q,u=transfer(name(ONE))
        s,v=transfer(name(B))
        I,W=[PchipInterpolator(np.log(x),y[:,0]) for x,y in ((q,u),(s,v))]
        grid=np.geomspace(.015,10.,501)
        probes=np.array([4.,8.,10.,16.,20.,22.])
        assert s[-1]>=22
        dif=W(np.log(grid))-I(np.log(grid))
        j=int(np.argmax(np.abs(dif)))
        ind=np.flatnonzero((s<22)&(v[:,0]<.5))
        if len(ind):
            m=int(ind[0])
            assert m and v[m-1,0]>.5
            half=float(np.exp(brentq(lambda r:float(W(r))-.5,np.log(s[m-1]),np.log(s[m]))))
        else:
            half=None
        rows.append(dict(z=z,probes=probes.tolist(),wide_higher_delta_over_cdm=W(np.log(probes)).tolist(),
            max_overlap_abs_defect=float(abs(dif[j])),overlap_k_at_max_h_per_Mpc=float(grid[j]),
            absolute_difference_at_k10=float(dif[-1]),wide_fluid_first_half_h_per_Mpc=half))
    return dict(changes=D,data=rows,scope="Closure retained; later switch, higher lmax and wider native k grid changed together. Any half crossing remains numerical and is not a physical transfer cutoff.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("native",type=Path)
    p.add_argument("--prepare",type=Path)
    p.add_argument("--cases",choices=list(LABELS),nargs="+")
    p.add_argument("--time-probe",action="store_true")
    p.add_argument("--wide",action="store_true")
    a=p.parse_args()
    if a.prepare is not None:
        result=prepare(a.native.resolve(),a.prepare.resolve())
    elif a.time_probe:
        result=history(a.native.resolve())
    elif a.wide:
        result=wider_default(a.native.resolve())
    else:
        result=compare(a.native.resolve(),a.cases)
    print(json.dumps(result,indent=2,allow_nan=False))
