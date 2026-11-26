"""Native one-versus-two species test of the same total phase-space input."""
import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator

from archived_collisionless_streaming_reference import Background
from class_exact_coordinate_remap_reference import scan
from conditional_relic_normalization_reference import conversion
from export_finite_cooling_class_psd_reference import smooth
from export_finite_cooling_pair_mixture_reference import support
from three_species_decay_onsager_reference import mp

PREFIX="finite_cooling_pairmix_tracer_ratio5_s14"
ONE=PREFIX+"_sync"
TWO=PREFIX+"_split_sync"
PROMPT="finite_cooling_tracer_ratio5_s14_sync"


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def fields(path):
    return {k.strip():v.strip() for k,v in (l.split("=",1) for l in path.read_text().splitlines()
                                          if l.strip() and not l.lstrip().startswith("#"))}


def prepare(eos,pair,joined,prompt,dest):
    raw=json.loads(pair.read_text())
    assert sha(eos)==raw["eos_sha256"]
    bg=Background(eos)
    qfac=bg.scale(.001)*.001/bg.T0
    fac=(2*np.pi**2/45)*bg.gs0/(mp*conversion(eos)["Omega_h2_per_GeV_times_yield"])
    v=raw["episodes"][-1]["modes"][-1]
    M=v["new_scalar_recoil_shape"]
    omega=float(np.sum(np.asarray(v["new_scalar_omega"])[:,0]))
    p,w=support(np.sum(M["count_fraction_by_epoch"],axis=0),np.sum(M["q2_weighted_count_fraction_by_epoch"],axis=0),fac*omega)
    q,f=smooth(qfac*p,w,.14,5200)
    f*=w.sum()/(4*np.pi*np.trapezoid(q*q*f,q))
    dest.mkdir(parents=True,exist_ok=True)
    fs=dest/"finite_pair_only_ratio5_s14.dat"
    np.savetxt(fs,np.c_[q,f],fmt="%.16e")
    for n in (160,224,320):
        control=scan(fs,n)["best"]
        if control["maximum_residual"]<1e-7:
            break
    assert control["maximum_residual"]<1e-7
    S=control["scale"]
    hot=dest/"finite_pair_only_ratio5_s14_remapped.dat"
    np.savetxt(hot,np.c_[q/S,f],fmt="%.16e")
    P=fields(prompt/(PROMPT+".ini"))
    J=fields(joined/(ONE+".ini"))
    match=("output","P_k_max_h/Mpc","z_pk","h","T_cmb","omega_b","A_s","n_s","tau_reio","N_ur",
           "gauge","omega_cdm","l_max_ncdm","ncdm_fluid_trigger_tau_over_tau_k")
    assert all(P[k]==J[k] for k in match)
    assert abs((float(P["omega_ncdm"])+omega)/float(J["omega_ncdm"])-1.)<1e-13
    P["root"]=str(dest/(TWO+"_"))
    P["N_ncdm"]="2"
    extra={"m_ncdm":86,"T_ncdm":S,"deg_ncdm":1,"use_ncdm_psd_files":1,
           "ncdm_psd_filenames":hot,"omega_ncdm":omega,"ncdm_quadrature_strategy":2,
           "ncdm_maximum_q":20,"ncdm_N_momentum_bins":n}
    for k,v in extra.items():
        P[k]=f"{P[k]},{v}"
    (dest/(TWO+".ini")).write_text("".join(f"{k} = {v}\n" for k,v in P.items()))
    return dict(pair_input_record_sha256=sha(pair),raw_hot_psd_sha256=sha(fs),mapped_hot_psd_sha256=sha(hot),
                hot_Omega_h2=omega,hot_momentum_nodes=n,external_hot_grid_screen=control,
                scope="Same scalar f0 split into two equal-mass tracers, not new physics. Native settings and total fixed density match the joined test.")


def indices(p,match):
    with p.open() as f:
        l=next(x for x in f if x.startswith("#") and match in x)
    return {k:int(i)-1 for i,k in re.findall(r"(\d+):(\S+)",l)}


def transfer(p):
    k=indices(p,"d_ncdm[0]")
    x=np.loadtxt(p)
    assert np.all(x[:,k["d_cdm"]]<0.)
    return x[:,0],np.column_stack([x[:,i]/x[:,k["d_cdm"]] for key,i in k.items() if key.startswith("d_ncdm[")])


def audit(joined,split,prompt,joined_psd,prompt_psd):
    extra=split/"finite_pair_only_ratio5_s14.dat"
    A=[np.loadtxt(f).T for f in (joined_psd,prompt_psd,extra)]
    x=np.unique(np.concatenate([a[0] for a in A]))
    d=np.array([np.nan_to_num(CubicSpline(q,F,extrapolate=False)(x),nan=0.) for q,F in A])
    err=[float(np.trapezoid(x**(q+2)*abs(d[0]-d[1]-d[2]),x)/np.trapezoid(x**(q+2)*d[0],x)) for q in (0,1,2,4)]
    assert max(err)<2e-7
    file=split/(TWO+"_00_background.dat")
    k=indices(file,"rho_ncdm[0]")
    a=np.loadtxt(file,usecols=[k[f] for f in ("z","(.)rho_ncdm[0]","(.)rho_ncdm[1]")])
    frac=a[:,2]/(a[:,1]+a[:,2])
    rows=[]
    for i,z in enumerate((0,5,20),1):
        w=float(np.interp(np.log1p(z),np.log1p(a[::-1,0]),frac[::-1]))
        paths=[r/f"{t}_00_z{i}_tk.dat" for r,t in ((joined,ONE),(split,TWO),(prompt,PROMPT))]
        cur=[transfer(p) for p in paths]
        kk=np.unique(np.r_[np.geomspace(.015,10,420),[.3,1.,2.,4.,8.,10.]])
        ss=[PchipInterpolator(np.log(q),S,axis=0)(np.log(kk)) for q,S in cur]
        s=ss[1]@np.array([1-w,w])
        e=ss[0][:,0]-s
        j=int(np.argmax(abs(e)))
        rows.append(dict(z=z,hot_energy_fraction=w,prompt_native_k10=float(ss[2][-1,0]),
            joined_k10=float(ss[0][-1,0]),cold_part_split_k10=float(ss[1][-1,0]),
            hot_part_split_k10=float(ss[1][-1,1]),energy_weighted_split_k10=float(s[-1]),
            joined_minus_split_k10=float(e[-1]),max_abs_joined_minus_split=float(abs(e[j])),
            k_at_max_h_per_Mpc=float(kk[j]),max_prompt_split_vs_separate=float(np.max(abs(ss[1][:,0]-ss[2][:,0])))))
    return dict(dense_input_moment_weighted_L1_errors=err,moment_orders=[0,1,2,4],rows=rows,
        scope="Metamorphic test of Boltzmann tracer numerical representation. Split density contrasts weighted by actual native scalar component energy. Independent collisionless hierarchies at fixed metric should add; fluid closure calibrated to a joined shape need not. Neither approximate closure certifies all-scale transfer.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for k in ("joined","split","prompt"):
        p.add_argument(k,type=Path)
    p.add_argument("--eos",type=Path)
    p.add_argument("--pair",type=Path)
    p.add_argument("--joined-psd",type=Path)
    p.add_argument("--prompt-psd",type=Path)
    a=p.parse_args()
    if a.joined_psd is not None:
        result=audit(a.joined,a.split,a.prompt,a.joined_psd,a.prompt_psd)
    else:
        result=prepare(a.eos,a.pair,a.joined,a.prompt,a.split)
    print(json.dumps(result,indent=2,allow_nan=False))
