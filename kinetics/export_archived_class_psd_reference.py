"""Numerical CLASS PSD bridge for finite-shell, post-production archives."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.ndimage import gaussian_filter1d

from archived_scalar_pair_reaction_reference import points,snapshot_functions
from archived_collisionless_streaming_reference import Background
from conditional_relic_normalization_reference import conversion
from three_species_decay_onsager_reference import m,mp


def export(p,w,ai,T0,mass,sigma,nx,out):
    # Physical number measure today, with actual spin counts already in w.
    q=p*ai/T0
    prob=w/w.sum()
    c=np.cumsum(prob)
    ii=(c>=1e-13)&(np.r_[0.,c[:-1]]<=1.-1e-13)&(prob>0.)
    assert np.sum(prob[~ii])<3e-12
    x=np.linspace(np.log(q[ii].min())-.55,np.log(q[ii].max())+.55,nx)
    dx=x[1]-x[0]
    u=(np.log(q[ii])-x[0])/dx
    j=np.floor(u).astype(int)
    v=u-j
    mom=np.bincount(j,prob[ii]*(1-v),minlength=nx)+np.bincount(j+1,prob[ii]*v,minlength=nx)
    dens=gaussian_filter1d(mom,sigma/dx,mode="constant",truncate=12.)
    qp=np.exp(x)
    F=dens/(qp**3*dx)
    # Make endpoint derivatives finite for CLASS's positive exponential
    # tail; floor density is numerically negligible after normalization.
    F+=F.max()*1e-18*np.exp(-qp/8)
    F*=w.sum()*(ai/T0)**3/(4*np.pi)
    F*=w.sum()*(ai/T0)**3/(4*np.pi*np.trapezoid(qp**2*F,qp))
    cs=CubicSpline(qp,F)
    tst=(qp[:-1,None]*np.array([.8,.5,.2])+qp[1:,None]*np.array([.2,.5,.8])).ravel()
    assert np.min(cs(tst))>0. and F[-1]<F[-2]
    ratio=np.array([np.trapezoid(qp**(2+k)*F,qp)/np.trapezoid(qp**2*F,qp)/(prob@q**k)
                    for k in (1,2,3)])
    assert max(abs(ratio-1.))<.0012
    np.savetxt(out,np.column_stack([qp,F]),fmt="%.16e")
    return dict(file=out.name,mass_eV=mass*1e9,q_max_recommended=float(qp[-1]),
                rest_number_per_today_Tgamma_cube=float(w.sum()*(ai/T0)**3),
                table_q_moments_1_2_3_over_shell=ratio.tolist(),clipped_probability=float(prob[~ii].sum()))


def generate(s,tab,out,sigma,nx):
    bg=Background(tab)
    ai=float(bg.scale(s["T_visible_GeV"]))
    p=np.asarray(s["lab_momentum_edges_GeV"])
    fp,_,fc,_=snapshot_functions(s)
    g0=2*np.pi**2/45*bg.gs0
    fac=conversion(tab)["Omega_h2_per_GeV_times_yield"]
    result=[]
    for name,mi,gi,f in (("stable",m,4,fc),("scalar",mp,1,fp)):
        xx,_,nn=points(p,mi,f,40)
        w=gi*nn
        row=export(xx,w,ai,bg.T0,mi,sigma,nx,out/(name+".dat"))
        row["omega_h2_from_retained_rest_count"]=row["rest_number_per_today_Tgamma_cube"]*mi*fac/g0
        result.append(row)
    return dict(snapshot_a=ai,sigma_log_q=sigma,grid_nodes=nx,components=result,
                convention="q=p_comoving/Tgamma_today, T_ncdm=1, deg_ncdm=1. File f0 includes spin counts divided by (2 pi)^3; exterior vacuum. Positive log smoothing exactly normalizes number and audits moment drift.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("table");p.add_argument("snapshot");p.add_argument("out")
    p.add_argument("--sigma",type=float,default=.010)
    p.add_argument("--nodes",type=int,default=3200)
    a=p.parse_args()
    out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    d=json.loads(Path(a.snapshot).read_text())["snapshots"][-1]
    r=generate(d,a.table,out,a.sigma,a.nodes)
    (out/"manifest.json").write_text(json.dumps(r,indent=2,allow_nan=False)+"\n")
    print(json.dumps(r,indent=2,allow_nan=False))
