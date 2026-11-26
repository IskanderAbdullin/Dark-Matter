"""First-order inverse-removal sensitivity on the same fixed metric."""
import argparse
import hashlib
import json
from pathlib import Path
import tempfile

import numpy as np
from scipy.interpolate import CubicSpline

from archived_collisionless_streaming_reference import Background
from conditional_relic_normalization_reference import conversion
from export_finite_cooling_class_psd_reference import smooth
from export_finite_cooling_pair_mixture_reference import support
from finite_cooling_relativistic_los_reference import run as los
from three_species_decay_onsager_reference import mp


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(eos,absorption,joined,original_los,metrics):
    data=json.loads(absorption.read_text())
    original=json.loads(original_los.read_text())
    assert data["eos_sha256"]==original["eos_sha256"]==sha(eos)
    assert original["spectrum_sha256"][0]==sha(joined)
    bg=Background(eos)
    factor=(2*np.pi**2/45)*bg.gs0/(mp*conversion(eos)["Omega_h2_per_GeV_times_yield"])
    B=float(bg.scale(.001))*.001/bg.T0
    am=mp/(21*bg.T0)
    q0,f0=np.loadtxt(joined).T
    W0=np.trapezoid(f0*q0*q0*np.hypot(q0,am),q0)
    reg=CubicSpline(q0,f0,extrapolate=False)
    report=[]
    cdf=[]
    for a in data["cases"]:
        h=a["removed_prompt_shape"]
        yy=np.asarray(h["absorbed_omega_by_epoch_and_bin"])
        zz=np.asarray(h["absorbed_q2_times_omega_by_epoch_and_bin"])
        U=float(yy.sum())
        p,m=support(yy.sum(axis=0)/U,zz.sum(axis=0)/U,U*factor)
        q,f=smooth(B*p,m,.14,5200)
        f*=m.sum()/(4*np.pi*np.trapezoid(f*q*q,q))
        assert abs(np.trapezoid(f*q**4,q)/(np.trapezoid(f*q*q,q)*B*B)-zz.sum()/U)<1e-8
        w=np.trapezoid(f*q*q*np.hypot(q,am),q)
        assert 0<w/W0<.003
        I=CubicSpline(q,f,extrapolate=False)
        grid=np.unique(np.r_[q,q0])
        diff=np.nan_to_num(reg(grid),nan=0.)-np.nan_to_num(I(grid),nan=0.)
        neg=np.trapezoid(np.maximum(0.,-diff)*grid**2,grid)/np.trapezoid(f0*q0*q0,q0)
        assert neg<2e-12
        with tempfile.TemporaryDirectory(prefix="scalar-removal-") as d:
            fp=Path(d)/"absorption.dat"
            np.savetxt(fp,np.c_[q,f],fmt="%.16e")
            rem=los(eos,[fp,fp,fp],metrics)
        rows=[]
        for O,R in zip(original["results"],rem["results"]):
            assert O["metric_sha256"]==R["metric_sha256"]
            k=O["k_h_per_Mpc"]
            key="los_over_cdm_joined_prompt_recoil_losIBP_direct_monopoleIC"
            base=O["cases"][-1][key][0][0]
            sub=R["cases"][-1][key][0][0]
            corrected=(base*W0-sub*w)/(W0-w)
            rows.append(dict(k_h_per_Mpc=k,conditional_absorbed_population_delta_over_cdm=sub,
                previous_forward_only_tracer_delta_over_cdm=base,after_passive_subtraction_delta_over_cdm=corrected,
                signed_absolute_change=corrected-base))
        cdf.append(np.cumsum(yy.sum(axis=0))/U)
        report.append(dict(control=a["configuration"],absorption_rest_omega_h2=U,
            removed_mixed_energy_fraction_at_z20=float(w/W0),worst_tested_negative_number_fraction=float(neg),z20=rows))
    return dict(input_sha256=[sha(p) for p in (eos,absorption,joined,original_los,*metrics)],
        max_resolved_absorption_CDF_between_controls=float(np.max(abs(cdf[0]-cdf[-1]))) if len(cdf)>1 else None,
        scope="Signed passive subtraction on the same fixed metric; inverse removes stocked prompt particles without recomputing the forward source or evolving feedback. Positivity checked on archive knot union, not rigorous positivity between knots. Not a coupled Boltzmann or total matter result.",cases=report)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for k in ("eos","absorption","joined","original_los"):
        p.add_argument(k,type=Path)
    p.add_argument("metrics",type=Path,nargs="+")
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.absorption,a.joined,a.original_los,a.metrics),indent=2,allow_nan=False))
