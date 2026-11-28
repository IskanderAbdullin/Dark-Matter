"""Positive first-order prompt plus spontaneous-pair reconstruction."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

from archived_collisionless_streaming_reference import Background
from conditional_relic_normalization_reference import conversion
from export_finite_cooling_class_psd_reference import smooth
from finite_cooling_reheat_portal_reference import TF
from three_species_decay_onsager_reference import mp


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def support(n,p2,target):
    n,p2=np.asarray(n),np.asarray(p2)
    assert np.all(n>=0.) and np.all(p2>=0.)
    k=(n>1e-15)&(p2>0.)
    assert n[~k].sum()<1e-12 and p2[~k].sum()<1e-11*p2.sum()
    return np.sqrt(p2[k]/n[k]),target*n[k]/n[k].sum()


def run(eos,prompt,pairs,out,nodes):
    a,b=[json.loads(f.read_text()) for f in (prompt,pairs)]
    esh=sha(eos)
    assert a["eos_sha256"]==b["eos_sha256"]==esh
    assert b["episodes"][0]["modes"][-1]["new_scalar_recoil_shape"] is not None
    bg=Background(eos)
    B=float(bg.scale(TF))*TF/bg.T0
    U=(2*np.pi**2/45)*bg.gs0/(mp*conversion(eos)["Omega_h2_per_GeV_times_yield"])
    out.mkdir(parents=True,exist_ok=True)
    summary=[]
    for p,s in zip(a["cases"],b["episodes"]):
        assert abs(p["T_R_GeV"]/s["T_R_GeV"]-1.)<1e-13
        v=s["modes"][-1]
        sm=v["new_scalar_recoil_shape"]
        pg=np.sum(p["count_fraction_per_epoch_species_bin"],axis=0)[0]
        p2=np.sum(p["p_over_T_squared_fractional_number_moment_per_epoch_species_bin"],axis=0)[0]
        S=np.sum(sm["count_fraction_by_epoch"],axis=0)
        S2=np.sum(sm["q2_weighted_count_fraction_by_epoch"],axis=0)
        n0=p["species"][0]["Omega_h2"]
        n1=np.sum(np.array(v["new_scalar_omega"])[:,0])
        q0,g0=support(pg,p2,U*n0)
        q1,g1=support(S,S2,U*n1)
        q=np.r_[q0,q1]*B
        n=np.r_[g0,g1]
        source_rms=np.sqrt((n@(q*q))/n.sum())/B
        assert abs(source_rms/v["first_order_mixture_scalar_rms_p_over_TF"]-1.)<2e-5
        shape="instant" if p["T_R_GeV"]==p["Tmax_GeV"] else "ratio5"
        for sigma in (.10,.14):
            xx,f=smooth(q,n,sigma,nodes)
            r=n.sum()/(4*np.pi*np.trapezoid(xx*xx*f,xx))
            f*=r
            m2=4*np.pi*np.trapezoid(xx**4*f,xx)
            assert abs(m2/(n@(q*q))-1.)<1e-7
            a1=(xx[:-1,None]*np.array([.8,.5,.2])+xx[1:,None]*np.array([.2,.5,.8])).ravel()
            assert np.min(CubicSpline(xx,f)(a1))>0 and 0<f[-1]<f[-2]
            fp=out/f"{shape}-scalar-pair-mix-sigma{sigma:.2f}.dat"
            np.savetxt(fp,np.column_stack([xx,f]),fmt="%.16e")
            summary.append(dict(case=shape,species="scalar",file=fp.name,sha256=sha(fp),
                actual_fixed_count_Omega_h2=n0+n1,prompt_omega_h2=n0,first_order_pair_omega_h2=n1,
                sigma_ln_q=sigma,reconstruction_anchor_p_over_T_rms=float(source_rms),
                relative_second_moment_defect=float(m2/(n@(q*q))-1.)))
    return dict(eos_sha256=esh,prompt_sha256=sha(prompt),pairs_sha256=sha(pairs),q_from_p_over_T=B,
       components=summary,scope="Positive provisional first-order shape at 1 MeV: prompt sample plus spontaneous pairs before the anchor, no depleted stable source or feedback, no pair reactions after anchor. Normalizes its own counts, not a coupled physical outcome.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for key in ("eos","prompt","pairs","out"):
        p.add_argument(key,type=Path)
    p.add_argument("--nodes",type=int,default=5200)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.prompt,a.pairs,a.out,a.nodes),indent=2,allow_nan=False))
