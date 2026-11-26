"""Controlled, positive smoothing of the recorded finite-cooling PSD moments."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import logsumexp

from archived_collisionless_streaming_reference import Background
from conditional_relic_normalization_reference import conversion
from finite_cooling_reheat_portal_reference import TF
from three_species_decay_onsager_reference import m, mp


def checksum(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def smooth(p,n,sigma,grid):
    """Lognormal with median correction enforcing the exact global p²."""
    logp=np.log(p)
    xc=logp-sigma**2
    u=np.linspace(xc.min()-12*sigma,xc.max()+12*sigma,grid)
    lp=logsumexp(np.log(n)-.5*((u[:,None]-xc)/sigma)**2,axis=1)
    f=np.exp(lp-3*u)/(4*np.pi*np.sqrt(2*np.pi)*sigma)
    return np.exp(u),f


def write(eos,distribution,out,nodes):
    raw=json.loads(Path(distribution).read_text())
    assert raw["eos_sha256"]==checksum(eos)
    b=Background(eos)
    factor=(2*np.pi**2/45)*b.gs0/conversion(eos)["Omega_h2_per_GeV_times_yield"]
    q_factor=float(b.scale(TF))*TF/b.T0
    root=Path(out)
    root.mkdir(parents=True,exist_ok=True)
    result=[]
    for case in raw["cases"]:
        val="instant" if case["Tmax_GeV"]==case["T_R_GeV"] else "ratio5"
        N=np.sum(case["count_fraction_per_epoch_species_bin"],axis=0)
        P2=np.sum(case["p_over_T_squared_fractional_number_moment_per_epoch_species_bin"],axis=0)
        for j,(mass,name) in enumerate(((mp,"scalar"),(m,"stable"))):
            allprob,allp2=N[j],P2[j]*q_factor**2
            keep=(allprob>1e-15)&(allp2>0.)
            assert np.sum(allprob[~keep])<1e-12 and np.sum(allp2[~keep])/allp2.sum()<1e-11
            p=np.sqrt(allp2[keep]/allprob[keep])
            omega=case["species"][j]["Omega_h2"]
            n=allprob[keep]*(factor*omega/(mass*allprob[keep].sum()))
            expect=[float(n.sum()),float(n@(p*p))]
            for sigma in (.10,.14):
                q,f=smooth(p,n,sigma,nodes)
                fr=n.sum()/(4*np.pi*np.trapezoid(q*q*f,q))
                f*=fr
                cs=CubicSpline(q,f)
                targets=(q[:-1,None]*np.array([.8,.5,.2])+q[1:,None]*np.array([.2,.5,.8])).ravel()
                assert np.min(cs(targets))>0 and 0<f[-1]<f[-2]
                evald=[float(4*np.pi*np.trapezoid(q**i*f,q)) for i in (2,4)]
                defect=np.array(evald)/expect-1.
                assert np.max(np.abs(defect))<1e-7
                fname=f"{val}-{name}-sigma{sigma:.2f}.dat"
                path=root/fname
                np.savetxt(path,np.column_stack([q,f]),fmt="%.16e")
                result.append(dict(case=val,species=name,mass_eV=mass*1e9,actual_fixed_count_Omega_h2=omega,
                     file=fname,sha256=checksum(path),sigma_ln_q=sigma,
                     expected_number_and_q2_moment_per_Tgamma3=expect,
                     relative_table_trapezoid_defects=defect.tolist(),
                     support_q_min_max=q[[0,-1]].tolist(),clipped_number_fraction=float(np.sum(allprob[~keep])),
                     number_weighted_q4_smoothing_factor=float(np.exp(4*sigma*sigma))))
    return dict(eos_sha256=checksum(eos),distribution_sha256=checksum(distribution),q_from_p_over_T=q_factor,
                components=result,scope="Numerical positive moment reconstruction of the recorded bins. It preserves total number and squared momentum, and is not an exact within-bin shape, a transfer computation, or a coupled production history.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos")
    p.add_argument("distribution")
    p.add_argument("out")
    p.add_argument("--nodes",type=int,default=5200)
    a=p.parse_args()
    print(json.dumps(write(a.eos,a.distribution,a.out,a.nodes),indent=2,allow_nan=False))
