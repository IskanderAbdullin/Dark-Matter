"""Independent additive late tracer Green function with a recorded metric drive."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_simpson
from scipy.interpolate import CubicSpline,PchipInterpolator

from archived_collisionless_streaming_reference import Background
from finite_cooling_pair_split_tracer_reference import indices
from three_species_decay_onsager_reference import mp


def run(eos,metric,files):
    bg=Background(eos)
    lab=indices(metric,"delta_ncdm[0]")
    tab=np.loadtxt(metric,usecols=[lab[k] for k in ("tau","a","psi")])
    k0=6.736
    with metric.open() as f:
        assert abs(float(f.readline().split("k = ")[1].split()[0])/k0-1)<1e-13
    _,j=np.unique(tab[:,1],return_index=True)
    tab=tab[j]
    ua=np.log(tab[:,1])
    ts=PchipInterpolator(ua,tab[:,0])
    vel=ts.derivative()
    psi=PchipInterpolator(ua,tab[:,2])
    raw=[np.loadtxt(p).T for p in files]
    N0=np.array([np.trapezoid(F*x*x,x) for x,F in raw])
    spline=[CubicSpline(q,F,extrapolate=False) for q,F in raw]
    floor=min(q[0] for q,_ in raw)
    ceiling=max(q[-1] for q,_ in raw)
    uf=np.log(1/21)
    out=[]
    for a0 in (.001,.005,.01):
        cut=min(ceiling,.10*mp*a0/bg.T0)
        q=np.geomspace(floor,cut,2300)
        w=np.array([np.maximum(0.,np.nan_to_num(s(q),nan=0.))*q*q for s in spline])
        N=np.trapezoid(w,q,axis=1)
        qtail=np.geomspace(cut,ceiling,1500)
        loss=np.array([np.trapezoid(np.maximum(0.,np.nan_to_num(s(qtail),nan=0.))*qtail*qtail,qtail) for s in spline])/N0
        assert np.max(loss)<1e-7
        S=[]
        for n in (380,970):
            u=np.linspace(np.log(a0),uf,n)
            a=np.exp(u)
            dv=vel(u)
            x=psi(u)
            assert np.all(dv>0.) and np.all(x*x[0]>0.)
            # tau derivative dz/dln(a), then ds=d tau/a.
            sn=-cumulative_simpson((dv/a)[::-1],x=-u[::-1],initial=0.)[::-1]
            travel=-sn
            assert travel.max()>0. and travel.min()>=0.
            weight=a*dv*x*travel
            ker=np.sinc((k0*bg.T0/(np.pi*mp))*q[:,None]*travel[None,:])
            norm=np.trapezoid(weight,u)
            resp=np.trapezoid(ker*weight,u,axis=1)/norm
            d=np.trapezoid(w*resp[None,:],q,axis=1)/N
            weighted=np.average(d[1:],weights=N[1:])
            es=abs(d[0]-weighted)
            assert es<4e-7 and np.min(d)>-.218 and np.max(d)<1.0000001
            S.append((d,weighted,es,float(norm)))
        residual=max(abs(S[0][0]-S[1][0]))
        assert residual<4e-5
        v,p,e,n0=S[1]
        out.append(dict(a_start=a0,z_final=20,k_Mpc_minus1=k0,
            normalized_joined_prompt_recoil_late_response=v.tolist(),weighted_split_response=float(p),
            joined_split_difference=float(e),time_quadrature_defect=float(residual),
            removed_number_fraction=loss.tolist(),
            positive_or_negative_metric_drive_sign=float(np.sign(psi(np.log(a0))))))
    return dict(cases=out,input_sha256=[hashlib.sha256(p.read_bytes()).hexdigest() for p in (eos,metric,*files)],
       scope="Newtonian nonrelativistic Green function, zero tracer perturbation at chosen a_start, same recorded metric potential drives every synthetic component. Uses no feedback or earlier hierarchy data, not the native evolving tracer or the physical total-matter transfer.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for k in ("eos","metric","joined","prompt","recoil"):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.metric,[a.joined,a.prompt,a.recoil]),indent=2,allow_nan=False))
