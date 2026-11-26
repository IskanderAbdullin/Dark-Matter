"""Retained low-reheat 86-eV scalar as a subdominant CLASS tracer."""
import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator
from scipy.optimize import brentq
from scipy.ndimage import gaussian_filter1d

from archived_class_manual_momentum_gate_reference import audit


def generate(fine,out):
    q,f=np.loadtxt(fine/"scalar.dat").T
    u=np.log(q)
    assert np.ptp(np.diff(u))<1e-12 and np.min(f)>0
    sigma=.055
    g=q*np.exp(-sigma*sigma)
    w=gaussian_filter1d(q**3*f,sigma/(u[1]-u[0]),mode="constant",truncate=10)/g**3
    w+=w.max()*1e-18*np.exp(-g/24)
    w*=np.trapezoid(q*q*f,q)/np.trapezoid(g*g*w,g)
    moments=[float(np.trapezoid(g**(j+2)*w,g)/np.trapezoid(q**(j+2)*f,q)) for j in range(5)]
    assert abs(moments[2]-1)<6e-6
    assert np.min(CubicSpline(g,w)((g[1:]+g[:-1])/2))>0 and w[-1]<w[-2]
    psd=out/"scalar-balanced.dat"
    np.savetxt(psd,np.c_[g,w],fmt="%.16e")
    J=json.loads((fine/"manifest.json").read_text())["components"]
    om=float(J[0]["omega_h2_from_retained_rest_count"])+1e-8
    omega=float(J[1]["omega_h2_from_retained_rest_count"])
    assert omega/om<3e-9 and omega>0.
    cases=[
      ("scalar_tracer_smooth128_sync",128,10,12,psd,1,18,96.),
      ("scalar_tracer_smooth256_sync",256,14,20,psd,1,18,96.),
      ("scalar_tracer_raw384_sync",384,14,20,fine/"scalar.dat",1,18,96.),
      ("scalar_tracer_raw160_q638_sync",160,14,20,fine/"scalar.dat",1,18,63.8376753507014),
      ("scalar_tracer_smooth128_wide_sync",128,10,12,psd,12,14,96.),
      ("scalar_tracer_smooth128_wide_delayed_sync",128,14,20,psd,12,14,96.)]
    for name,n,l,trigger,path,kmax,nk,qmax in cases:
        s=f"""output = mPk,dTk
P_k_max_h/Mpc = {kmax}
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
k_per_decade_for_pk = {nk}
tol_perturbations_integration = 1e-07
root = {out}/{name}_
omega_cdm = {om-omega:.17g}
N_ncdm = 1
m_ncdm = 86
T_ncdm = 1
deg_ncdm = 1
use_ncdm_psd_files = 1
ncdm_psd_filenames = {path}
tol_ncdm = 1e-05
tol_ncdm_bg = 1e-07
omega_ncdm = {omega:.17g}
ncdm_quadrature_strategy = 3
ncdm_maximum_q = {qmax:.17g}
ncdm_N_momentum_bins = {n}
ncdm_fluid_trigger_tau_over_tau_k = {trigger}
l_max_ncdm = {l}
k_min_tau0 = 100
"""
        (out/(name+".ini")).write_text(s)
    return dict(scalar_density_omega_h2=omega,dark_density_omega_h2=om,
                scalar_fraction_of_dark_density=omega/om,smoothing_moment_ratios_upto4=moments,
                psd_sha256=[hashlib.sha256(p.read_bytes()).hexdigest() for p in (fine/"scalar.dat",psd)],
                smooth_manual_quadratures=audit(psd,[128,256],96.)["results"],
                raw_manual_quadratures=audit(fine/"scalar.dat",[128,384],96.)["results"],
                independent_shifted_raw160_moment_control=audit(fine/"scalar.dat",[160],63.8376753507014))


def density_ratio(path,z):
    with path.open() as f:
        h=f.readline()
        assert re.search(rf"redshift z={z}(?:[\s.]|$)",h)
        labels=next(x for x in f if x.startswith("#") and "d_ncdm[0]" in x and "d_cdm" in x)
    cols={name:int(idx)-1 for idx,name in re.findall(r'(\d+):(\S+)',labels)}
    a=np.loadtxt(path)
    den=a[:,cols["d_cdm"]]
    assert den.min()<0 and np.max(den)<0
    return a[:,0],a[:,cols["d_ncdm[0]"]]/den


def extract(out,cases,anchor=None):
    rows=[]
    for i,z in enumerate((0,5,20),1):
        grid=np.array([.015,.03,.1,.3,.7,1.])
        names=[]
        for x in cases:
            f=out/f"{x}_00_z{i}_tk.dat"
            q,r=density_ratio(f,z)
            assert np.min(q)<=grid.min()<=grid.max()<=np.max(q)
            interp=PchipInterpolator(np.log(q),r)
            values=interp(np.log(grid))
            extra=np.array([2.,4.,8.,10.])
            high=extra[extra<=q.max()]
            ix=np.flatnonzero((r[:-1]>.5)&(r[1:]<=.5))
            first=None if len(ix)==0 else float(np.exp(brentq(
                lambda logk: float(interp(logk))-.5,np.log(q[ix[0]]),np.log(q[ix[0]+1]))))
            names.append(dict(case=x,ratio_at_probes=values.tolist(),
                    first_sample_bracketed_amplitude_half_h_per_Mpc=first,
                    optional_high_k_h_per_Mpc=high.tolist(),
                    optional_high_k_ratio=interp(np.log(high)).tolist(),
                    native_k_h_per_Mpc=q.tolist(),native_delta_scalar_over_delta_cdm=r.tolist()))
        critical=np.array([a["ratio_at_probes"] for a in names])
        spread=np.ptp(critical,axis=0)
        rows.append(dict(z=z,k_h_per_Mpc=grid.tolist(),cases=names,
                         probe_max_variant_spread=float(max(spread))))
    details=[]
    if anchor is not None:
        for c in cases:
            p=out/f"{c}_00_background.dat"
            with p.open() as f:
                label=next(x for x in f if "rho_ncdm[0]" in x and "conf. time" in x)
            assert "3:conf. time [Mpc]" in label
            cols={name:int(idx)-1 for idx,name in re.findall(r'(\d+):(\S+)',label)}
            x=np.loadtxt(p,usecols=(0,2,3,cols["(.)rho_ncdm[0]"]))
            a=-np.log1p(x[:,0])
            tau=float(np.interp(np.log(anchor),a,x[:,1]))
            l=(out/(c+".log")).read_text()
            u=np.array(re.findall(r'Mode k=([.e+\-\d]+): initializing vector at tau=([.e+\-\d]+)',l),float)
            assert len(u)>20 and np.min(u[:,1])>tau
            w=.6736**2*x[-1,3]/x[-1,2]**2
            details.append(dict(case=c,tau_at_archive_Mpc=tau,
                  largest_mode_k_times_tau_at_archive=float(np.max(u[:,0])*tau),
                  largest_k_times_tau_at_native_perturbation_start=float(np.max(u[:,0]*u[:,1])),
                  present_NATIVE_energy_density_omega_h2=float(w)))
    return dict(rows=rows,anchor_diagnostics=details,scope="Scalar contrast relative to CDM in the same synchronous adiabatic CLASS run, with dominant stable matter treated as cold. Physical scalar abundance fixed at its tiny retained value. Not a prediction in which the scalar makes all dark matter.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("out",type=Path)
    p.add_argument("--fine",type=Path)
    p.add_argument("--compare",nargs="+")
    a=p.parse_args()
    a.out=a.out.resolve()
    a.out.mkdir(parents=True,exist_ok=True)
    if a.compare:
        anchor=None if a.fine is None else json.loads((a.fine/"manifest.json").read_text())["snapshot_a"]
        r=extract(a.out,a.compare,anchor)
    else:
        r=generate(a.fine.resolve(),a.out)
    print(json.dumps(r,indent=2,allow_nan=False))
