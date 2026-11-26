"""Pressure-preserving log-q smoothing control for archived CLASS inputs."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator
from scipy.ndimage import gaussian_filter1d

from archived_class_manual_momentum_gate_reference import audit


def smoothing(path, output, sigma=.055):
    q, f = np.loadtxt(path).T
    x = np.log(q)
    d = np.diff(x)
    assert d.min() > 0 and np.ptp(d) < 1e-12 and np.min(f) > 0
    measure = q**3 * f
    norm = np.trapezoid(f*q*q,q)
    m = np.array([np.trapezoid(q**(2+j) * f, q) / norm for j in range(1,5)])
    # A Gaussian log shift multiplies <q^j> by exp(j^2*sigma^2/2).
    # The compensating shift preserves the q^2 pressure moment.
    g = q * np.exp(-sigma**2)
    h = gaussian_filter1d(measure, sigma/d[0], mode="constant", truncate=10.) / g**3
    h += h.max() * 1e-18 * np.exp(-g/8)
    h *= norm / np.trapezoid(h*g*g,g)
    values = CubicSpline(g,h)((g[:-1,None]*np.array([.2,.5,.8])
                              +g[1:,None]*np.array([.8,.5,.2])).ravel())
    assert np.min(values) > 0 and h[-1] < h[-2]
    moments = np.array([np.trapezoid(g**(2+j)*h,g)/norm for j in range(1,5)])
    ratio = moments/m
    assert abs(ratio[1]-1.) < 6e-6
    np.savetxt(output,np.c_[g,h],fmt="%.16e")
    spec = [64,96,128,160,256]
    # Homogeneous moments use q_i^2 times f0, epsilon=sqrt(q^2+(am/T0)^2).
    # This is a tabulated-shape check, not a Friedmann or production run.
    T0_eV=8.617333262e-5*2.72548
    a=np.geomspace(1e-13,1.,90)
    def ratios(nodes, distribution):
        am=a[:,None]*2e7/T0_eV
        eps=np.hypot(nodes,am)
        wt=nodes**2*distribution
        return np.column_stack([np.trapezoid(eps*wt,nodes,axis=1),
                                np.trapezoid(wt*nodes**2/(3*eps),nodes,axis=1)])
    comparison=ratios(g,h)/ratios(q,f)-1.
    background=[dict(scale_factor=float(a[j]),relative_drift=dict(density=float(v[0]),pressure=float(v[1])))
                for j,v in enumerate(comparison)]
    return dict(sigma_log_q=sigma,
                input_and_output_sha256=[hashlib.sha256(Path(z).read_bytes()).hexdigest()
                                         for z in (path,output)],
                raw_over_output_number=float(norm/np.trapezoid(g**2*h,g)),
                q_moments_after_over_before=ratio.tolist(),
                maximum_background_density_and_pressure_drift=np.max(np.abs(comparison),axis=0).tolist(),
                background_moment_drift=background,
                raw_manual=audit(path,spec,20.)["results"],
                smoothed_manual=audit(output,spec,20.)["results"])


def config(out,name,dark_density,psd=None,n=None,high_sampling=False,trigger=8):
    anchor=1e-8
    conf=f"""output = mPk
P_k_max_h/Mpc = 3
z_pk = 0,5
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
perturbations_verbose = {3 if psd else 1}
write background = yes
headers = yes
k_per_decade_for_pk = {45 if high_sampling else 18}
tol_perturbations_integration = 1e-7
root = {out}/{name}_
"""
    if not high_sampling:
        conf += "k_min_tau0 = 100\n"
    if psd is None:
        conf += f"omega_cdm = {dark_density+anchor}\n"
    else:
        conf+=f"""omega_cdm = {anchor}
N_ncdm = 1
m_ncdm = 20000000
T_ncdm = 1
deg_ncdm = 1
use_ncdm_psd_files = 1
ncdm_psd_filenames = {psd}
tol_ncdm = 1e-5
tol_ncdm_bg = 1e-7
omega_ncdm = {dark_density:.17g}
ncdm_quadrature_strategy = 3
ncdm_maximum_q = 20
ncdm_N_momentum_bins = {n}
ncdm_fluid_trigger_tau_over_tau_k = {trigger}
l_max_ncdm = 8
"""
    (out/(name+".ini")).write_text(conf)


def read(out, case, redshift_index):
    table = np.loadtxt(out/f"{case}_00_z{redshift_index}_pk.dat")
    assert table.shape[1] == 2 and np.all(table > 0)
    assert np.all(np.diff(table[:,0]) > 0.)
    return table


def at(table,k):
    assert table[0,0] <= np.min(k) <= np.max(k) <= table[-1,0]
    return np.exp(PchipInterpolator(np.log(table[:,0]),np.log(table[:,1]))(np.log(k)))


def spectrum(out,others):
    rows=[]
    for i,z in enumerate((0,5),1):
        cd=read(out,"cdm_balanced_grid_sync",i)
        rr=read(out,"cdm_true_sync",i)
        b=read(out,"psd_balanced_sync",i)
        mask=(cd[:,0]>=.012)&(cd[:,0]<=2.9)
        k,P=cd[mask].T
        assert len(k)>20 and b.shape==cd.shape
        PB=at(b,k)
        delta0=np.max(np.abs(b[:,0]/cd[:,0]-1.))
        shift=np.max(np.abs(PB/b[mask,1]-1.))
        assert delta0<1e-6 and shift<1e-5
        trial={c:at(read(out,c,i),k)/P for c in others}
        rows.append(dict(z=z,mode_count=len(k),k_bounds_h_per_Mpc=[k[0],k[-1]],
                         max_relative_native_node_shift=float(delta0),
                         max_smooth_self_interpolation_effect=float(shift),
                         max_CDM_grid_interpolation_plus_solver_difference=float(np.max(np.abs(at(rr,k)/P-1.))),
                         archived_ratio_curve_on_reference_nodes=dict(k_h_per_Mpc=k.tolist(),
                             balanced_over_CDM=(PB/P).tolist(),
                             alternative_CDM_over_CDM=(at(rr,k)/P).tolist(),
                             other_cases_over_CDM={c:v.tolist() for c,v in trial.items()}),
                         max_pressure_balanced_P_over_matched_CDM_minus_one=float(np.max(np.abs(PB/P-1.))),
                         control_cases={c:dict(max_P_over_matched_CDM_minus_one=float(np.max(np.abs(v-1.))),
                                               max_pressure_balanced_minus_case_over_CDM=float(np.max(np.abs(PB/P-v))))
                                        for c,v in trial.items()}))
    return rows


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("out",type=Path)
    p.add_argument("--fine-export",type=Path)
    p.add_argument("--standard-export",type=Path)
    p.add_argument("--audit",type=Path)
    p.add_argument("--compare",nargs="+")
    a=p.parse_args()
    a.out=a.out.resolve()
    a.out.mkdir(parents=True,exist_ok=True)
    if a.compare is not None:
        print(json.dumps(dict(rows=spectrum(a.out,a.compare)),indent=2,allow_nan=False))
    else:
        raw=(a.fine_export/"stable.dat").resolve()
        ref=(a.standard_export/"stable.dat").resolve()
        def density(sub):
            return json.loads((sub/"manifest.json").read_text())["components"][0]["omega_h2_from_retained_rest_count"]
        om=density(a.fine_export)
        om0=density(a.standard_export)
        assert abs(om/om0-1.)<2e-6
        result=smoothing(raw,a.out/"balanced-0.055.dat")
        conf=(("cdm_balanced_grid_sync",om,None,None,False,8),
              ("cdm_true_sync",om,None,None,True,8),
              ("psd_balanced_sync",om,a.out/"balanced-0.055.dat",256,False,8),
              ("psd_true_sync_base",om,ref,64,True,4),
              ("psd_true_sync_fine",om,raw,128,True,8),
              ("psd_raw64_matched_sync",om,raw,64,False,8),
              ("psd_raw256_matched_sync",om,raw,256,False,8))
        for args in conf: config(a.out,*args)
        a.audit.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
