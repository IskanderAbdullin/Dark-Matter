"""Independent momentum-resolved collisionless LOS on recorded CLASS metrics."""
import argparse
import hashlib
import json
from pathlib import Path
import re

import numpy as np
from scipy.integrate import cumulative_simpson,simpson
from scipy.interpolate import CubicSpline,PchipInterpolator
from scipy.special import spherical_jn

from archived_collisionless_streaming_reference import Background
from finite_cooling_pair_split_tracer_reference import indices
from three_species_decay_onsager_reference import mp


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def trajectory(metric,redshift=20.,cdm_only=False):
    if metric.suffix==".npz":
        assert redshift==20. and not cdm_only
        with np.load(metric,allow_pickle=False) as z:
            k=float(z["k"])
            endpoints,columns,start,native=[z[q] for q in ("tau_endpoints","drive","initial","z20_native")]
            assert columns.ndim==2 and columns.shape[1]==3 and start.size==9 and native.size in (2,3)
            assert 300<len(columns)<30000 and np.all(np.isfinite(columns))
            assert np.all(np.diff(endpoints)>0.) and np.all(np.diff(columns[:,0])>0.)
            ts=np.geomspace(*endpoints,len(columns))
            return k,start,float(endpoints[1]),native,np.column_stack([ts,columns.astype(float)]),str(z["raw_sha256"])
    with metric.open() as f:
        k=float(re.search(r"k = (\S+)",f.readline())[1])
    ind=indices(metric,"delta_cdm" if cdm_only else "delta_ncdm[0]")
    keys=("tau","a","psi","phi","delta_ur","theta_ur","shear_ur","delta_cdm")
    if not cdm_only:
        keys+=("delta_ncdm[0]","cs2_ncdm[0]")
    M=np.loadtxt(metric,usecols=[ind[t] for t in keys])
    _,j=np.unique(M[:,0],return_index=True)
    M=M[j]
    start=M[0,:9]
    tf=float(PchipInterpolator(np.log(M[:,1]),M[:,0])(np.log(1/(1+redshift))))
    native=PchipInterpolator(M[:,0],M[:,7:])(tf)
    return k,start,tf,native,M[:,:4],sha(metric)


def pack(metrics,dest,n):
    dest.mkdir(parents=True,exist_ok=True)
    rows=[]
    assert 300<n<30000
    for m in metrics:
        k,start,tf,native,T,h=trajectory(m)
        assert m.suffix!=".npz"
        tau=np.geomspace(start[0],tf,n)
        cols=PchipInterpolator(T[:,0],T[:,1:],axis=0)(tau).astype(np.float32)
        p=dest/f"finite-cooling-native-metric-k{k:g}.npz"
        np.savez_compressed(p,k=k,tau_endpoints=np.array([start[0],tf]),initial=start,
            z20_native=native,drive=cols,raw_sha256=np.array(h))
        rows.append(dict(file=p.name,k_Mpc_minus1=k,bytes=p.stat().st_size,
            number_of_geometric_time_samples=n,raw_metric_sha256=h,fixture_sha256=sha(p)))
    return dict(files=rows,scope="Portable float32 background/metric-only fixture, not a native spectrum. Native initial data and endpoint reference values are preserved in float64.")


def run(eos,files,metrics,redshift=20.,cdm_only=False):
    bg=Background(eos)
    raw=[np.loadtxt(p).T for p in files]
    sp=[CubicSpline(q,f,extrapolate=False) for q,f in raw]
    lo=min(np.log(q[0]) for q,f in raw)
    hi=max(np.log(q[-1]) for q,f in raw)
    answer=[]
    for metric in metrics:
        k,A0,tf,native,M,origin=trajectory(metric,redshift,cdm_only)
        assert 0<k*A0[0]<.10
        theta=np.array([A0[4]/4.,A0[5]/(3*k),A0[6]/2])
        aG=PchipInterpolator(M[:,0],M[:,1:4],axis=0)
        dPhi=PchipInterpolator(M[:,0],M[:,3]).derivative()
        data=[]
        controls=((720,10000),(1300,23000)) if cdm_only else ((720,9000),(1300,19000))
        for nq,nt in controls:
            q=np.exp(np.linspace(lo,hi,nq))
            t=np.geomspace(A0[0],tf,nt)
            a,psi,phi=aG(t).T
            am=a*(mp/bg.T0)
            derivative=dPhi(t)
            finale=np.hypot(q,am[-1])
            values=[]
            # Chunk the real-space phase integrals; no three-dimensional
            # perturbation-history arrays or growing multipole hierarchy.
            for bb in range(0,nq,24):
                xx=q[bb:bb+24,None]
                v=xx/np.hypot(xx,am[None,:])
                travel=cumulative_simpson(v[:,::-1],x=-t[::-1],axis=1,initial=0.)[:,::-1]
                arg=k*travel
                j1=spherical_jn(1,arg)
                x0=arg[:,0]
                j0=spherical_jn(0,x0)
                source=simpson(j1*k*(v*phi+psi/v),x=t,axis=1)
                ini0=-theta[0]*j0
                ini=ini0+3*theta[1]*spherical_jn(1,x0)-5*theta[2]*spherical_jn(2,x0)
                # The monopole projection is the momentum derivative
                # times IC + gravitational response. IBP removes phi'.
                integrate_by_parts=ini-phi[-1]+phi[0]*j0+source
                direct=ini+simpson(-derivative*spherical_jn(0,arg)+k*(psi/v)*j1,x=t,axis=1)
                values.append(np.stack([integrate_by_parts,direct,integrate_by_parts-ini+ini0]))
            F=np.concatenate(values,axis=1)
            num=[];den=[];press=[];it0=[]
            eps0=np.hypot(q,am[0])
            for coeff in sp:
                w=np.nan_to_num(coeff(q),nan=0.)
                dq=q*np.nan_to_num(coeff(q,1),nan=0.)
                den.append(np.trapezoid(q*q*finale*w,q))
                num.append(np.trapezoid(q*q*finale*dq*F,q,axis=1))
                press.append(np.trapezoid(q**4/(3*finale)*dq*F,q,axis=1))
                it0.append(-theta[0]*np.trapezoid(q*q*eps0*dq,q)/np.trapezoid(q*q*eps0*w,q))
            D=np.array(num)/np.array(den)[:,None]/native[0]
            pres=np.array(press)/np.array(num)
            split=float(np.average(D[1:,0],weights=den[1:]))
            derivative_defect=float(max(abs(D[:,0]-D[:,1])))
            assert derivative_defect<(2e-7 if cdm_only else 2e-8),(float(k),nq,nt,derivative_defect)
            assert abs(D[0,0]-split)<3e-7
            data.append(dict(q_nodes=nq,tau_nodes=nt,
               los_over_cdm_joined_prompt_recoil_losIBP_direct_monopoleIC=D.tolist(),
               max_absolute_ibp_vs_direct_phi_dot=derivative_defect,
               summed_prompt_recoil_with_full_IC=split,native_raw_trace_delta_over_cdm=float(native[1]/native[0]) if native.size>=2 else None,
               independent_deltaP_over_deltaRho_joined_prompt_recoil=pres[:,0].tolist(),
               native_raw_deltaP_over_deltaRho=float(native[2]) if native.size==3 else None,
               initial_joined_LOS_vs_native_delta=float(it0[0]/A0[8]-1.) if A0.size>=9 else None,
               initial_joined_prompt_recoil_delta=it0))
        low,high=[np.array(d["los_over_cdm_joined_prompt_recoil_losIBP_direct_monopoleIC"]) for d in data]
        assert np.max(abs(low-high))<2e-6,(float(k),float(redshift),np.max(abs(low-high)))
        answer.append(dict(k_Mpc_minus1=k,k_h_per_Mpc=k/(bg.H0/100),tau_final_Mpc=tf,
             initial_k_tau=float(k*A0[0]),first_recorded_relative_density=float(A0[8]/A0[7]) if A0.size>=9 else None,
             metric_sha256=sha(metric),original_native_metric_sha256=origin,joined_kinematic_density_probe_relative_defect=float(high[0,0]/(native[1]/native[0])-1) if native.size>=2 else None,
             max_time_and_momentum_refinement_difference=float(np.max(abs(low-high))),cases=data))
    return dict(eos_sha256=sha(eos),spectrum_sha256=[sha(f) for f in files],z=int(redshift) if (not cdm_only and float(redshift).is_integer()) else redshift,results=answer,
        scope="Collisionless Newtonian-gauge line-of-sight density on a fixed recorded metric, with approximate early adiabatic monopole/dipole/quadrupole and identical histories for all scalar shapes. Independent time/angular propagation, not metric generation, nonlinear production or total matter transfer.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for k in ("eos","joined","prompt","recoil"):
        p.add_argument(k,type=Path)
    p.add_argument("metrics",nargs="+",type=Path)
    p.add_argument("--pack",type=Path)
    p.add_argument("--pack-nodes",type=int,default=2400)
    p.add_argument("--z",type=float,default=20.)
    p.add_argument("--cdm-only",action="store_true")
    a=p.parse_args()
    out=pack(a.metrics,a.pack,a.pack_nodes) if a.pack is not None else run(a.eos,[a.joined,a.prompt,a.recoil],a.metrics,a.z,a.cdm_only)
    print(json.dumps(out,indent=2,allow_nan=False))
