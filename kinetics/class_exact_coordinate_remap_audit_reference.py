"""Compare completed native CLASS runs under an exact PSD coordinate remap."""
import argparse
import json
from pathlib import Path
import re

import numpy as np
from scipy.interpolate import PchipInterpolator

PROBES=np.array([.09,.12,.18,.24,.30,.50,.60])


def background(out,name):
    p=out/f"{name}_00_background.dat"
    with p.open() as f:
        line=next(s for s in f if "rho_ncdm[0]" in s)
    ids={k:int(i)-1 for i,k in re.findall(r"(\d+):(\S+)",line)}
    ids=list(ids[k] for k in ("z","H","(.)rho_ncdm[0]","(.)p_ncdm[0]"))
    z,H,rho,press=np.loadtxt(p,usecols=ids).T
    a=1/(1+z)
    assert np.all(np.diff(a)>0) and min(press)>0 and min(rho)>0
    return .6736**2*rho[-1]/H[-1]**2,a,press/rho


def survey(out,names,stable_table=None):
    ref="cdm_short_altq_sync"
    background_rows=[]
    for c in names:
        om,a,w=background(out,c)
        aa=np.array([1e-11,1e-9,1e-7,1e-4,1.])
        assert a[0]<aa[0] and aa[-1]==a[-1]
        assert abs(om/.11816900374571636-1)<1e-10
        v=PchipInterpolator(np.log(a),np.log(w))
        rr=dict(case=c,actual_Omega_h2=om,scale_factors_for_pressure=aa.tolist(),
                pressure_over_energy=np.exp(v(np.log(aa))).tolist())
        if stable_table is not None:
            from archived_class_massive_moment_dilation_reference import spline_tail,native
            S=(out/(c+".ini")).read_text()
            match=lambda k:re.search(r"(?m)^"+re.escape(k)+r" = ([^\n]+)$",S)[1]
            assert int(match("ncdm_quadrature_strategy"))==2
            scale=float(match("T_ncdm"))
            n=int(match("ncdm_N_momentum_bins"))
            q,w=native(n,None)
            q,w=q*scale,w*scale
            _,f=spline_tail(stable_table)
            F,D=f(q)
            rr["scaled_native_momentum_nodes"]=n
            rr["coordinate_scale"]=scale
            rr["independent_spline_number_dilation_residual"]=float((w@(q**3*D))/(-3*(w@(q*q*F)))-1)
        background_rows.append(rr)
    result=[]
    for i,z in enumerate((0,5),1):
        def get(c):
            x=np.loadtxt(out/f"{c}_00_z{i}_pk.dat")
            assert np.all(x>0.)
            return x
        d=get(ref)
        files=[get(c) for c in names]
        lo=max(v[0,0] for v in [d,*files])
        hi=min(v[-1,0] for v in [d,*files])
        k=np.unique(np.r_[d[(d[:,0]>lo)&(d[:,0]<hi),0],PROBES])
        assert lo<min(k) and max(k)<hi
        def curve(t):
            return np.exp(PchipInterpolator(np.log(t[:,0]),np.log(t[:,1]))(np.log(k)))
        baseline=curve(d)
        vals=[]
        for v,c in zip(files,names):
            y=curve(v)/baseline
            m=int(np.argmax(abs(y-1.)))
            vals.append(dict(case=c,max_delta_P_over_CDM_minus_one=float(abs(y[m]-1.)),
                signed_delta_at_largest_probe=float(y[m]-1.),k_at_largest_probe_h_per_Mpc=float(k[m]),
                k_probe_h_per_Mpc=PROBES.tolist(),
                power_over_CDM_at_probes=y[np.searchsorted(k,PROBES)].tolist(),
                signed_common_grid_fractional_residual=(y-1.).tolist()))
        result.append(dict(z=z,comparison_k_grid_h_per_Mpc=k.tolist(),cases=vals))
    return dict(rows=result,background_checks=background_rows,scope="Native unmodified CLASS curves, same stable physical f0 but different numerical coordinate scale/number quadrature. Short wavenumber range. No perturbative production history or proof of physical small-scale transfer convergence.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("native",type=Path)
    p.add_argument("--stable-table",type=Path)
    p.add_argument("--cases",nargs="+",default=["psd_coordinate_native128_sync","psd_coordinate_remap128_sync"])
    a=p.parse_args()
    print(json.dumps(survey(a.native,a.cases,a.stable_table),indent=2,allow_nan=False))
