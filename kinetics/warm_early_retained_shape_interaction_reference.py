"""Matched warm retained shape sensitivity, with no guessed exterior spectrum."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from scalar_pair_expanding_continuation_reference import MASS, moments
from warm_early_serial_support_audit_reference import sha

PROBES=np.array([.003,.01,.03,.10,.3,1.,2.,3.,5.,10.])


def read(p):
    return json.loads(gzip.decompress(p.read_bytes()) if p.suffix==".gz" else p.read_bytes())


def integrate(p,th,mass,g,n,upper):
    a=p[:-1]
    b=p[1:]
    c=np.minimum(b,upper)
    j=np.flatnonzero(c>a)
    nodes,weights=roots_legendre(n)
    q=((a[j]+c[j])[:,None]+(c[j]-a[j])[:,None]*nodes)/2.
    dq=(c[j]-a[j])[:,None]*weights/2.
    E=np.hypot(q,mass)
    K=q*q/(E+mass)
    ke=p*p/(np.hypot(p,mass)+mass)
    t=th[0,j,None]+th[1,j,None]*(K-ke[j,None])/np.diff(ke)[j,None]
    if g==1:
        assert np.all(t>0.)
        f=np.exp(-t)/(-np.expm1(-t))
    else:
        from scipy.special import expit
        f=expit(-t)
    u=dq*q*q*(g/(2*np.pi**2))*f
    return np.array([u.sum(),(u*K).sum(),(u*q*q).sum(),(u*q**4).sum(),(u*q*q/(3*E)).sum()])


def properties(S):
    q=np.asarray(S["lab_momentum_edges_GeV"])
    th=[np.asarray(S[k]) for k in ("entropy_variable_scalar","entropy_variable_fermion")]
    T=S["T_visible_GeV"]
    target=moments(q,th).sum(axis=2)
    rows=[]
    for j,m in enumerate(MASS):
        ref=integrate(q,th[j],m,(1,4)[j],40,q[-1])
        fine=integrate(q,th[j],m,(1,4)[j],64,q[-1])
        d=float(np.max(abs(fine/ref-1.)))
        err=float(np.max(abs(fine[:2]/target[j]-1.)))
        assert max(d,err)<6e-8,(j,d,err)
        cdf=np.array([integrate(q,th[j],m,(1,4)[j],40,x*T)[0]/fine[0] for x in PROBES])
        assert min(cdf)>=0. and max(cdf)<=1.+5e-8 and np.min(np.diff(cdf))>=-1e-12
        rows.append(dict(number_GeV3=fine[0],mean_K_GeV=fine[1]/fine[0],rms_p_over_T=np.sqrt(fine[2]/fine[0])/T,
                         fourth_moment_p_over_T=fine[3]/fine[0]/T**4,mean_pressure_per_particle_GeV=fine[-1]/fine[0],
                         count_CDF=cdf.tolist(),relative_quadrature_and_moment_checks=[d,err]))
    return rows


def run(checkpoints,base=4):
    assert base in (4,8)
    if len(checkpoints)!=5:
        raise ValueError("Supply q16/4, q32/4, q16/8, q32/8, q16/32.")
    B=[read(p) for p in checkpoints]
    expected=[base,base,2*base,2*base,32]
    for n,z in zip(expected,B,strict=True):
        assert z["blocks_finished"]==z["blocks_total"]==n
        for key in ("schema","eos_sha256","local_prompt_jacobian","reference_sha256"):
            assert z["plan"][key]==B[0]["plan"][key]
    assert B[0]["plan"]["source_sha256"]==B[2]["plan"]["source_sha256"]==B[4]["plan"]["source_sha256"]
    assert B[1]["plan"]["source_sha256"]==B[3]["plan"]["source_sha256"]
    S=[x["snapshot"] for x in B]
    T=S[0]["T_visible_GeV"]
    assert max(abs(x["T_visible_GeV"]/T-1.) for x in S)<1e-12
    assert max(abs(x["scale_factor"]/S[0]["scale_factor"]-1.) for x in S)<1e-12
    data=[properties(x) for x in S]
    comparisons=[]
    for j,label in enumerate(("scalar","stable")):
        z=[x[j] for x in data]
        V=np.array([[x[k] for k in ("rms_p_over_T","fourth_moment_p_over_T","mean_pressure_per_particle_GeV")] for x in z])
        U=np.array([x["count_CDF"] for x in z])
        pair4=np.log(V[1]/V[0])
        pair8=np.log(V[3]/V[2])
        comparisons.append(dict(species=label,relative_moment_log_support_4_and_8=[pair4.tolist(),pair8.tolist()],
                                log_support_time_interaction=(pair8-pair4).tolist(),log_time_shift_q16_32_vs_8=np.log(V[4]/V[2]).tolist(),
                                maximum_CDF_support_difference_4_and_8=[float(max(abs(U[1]-U[0]))),float(max(abs(U[3]-U[2])))],
                                maximum_CDF_time_difference_q16_32_vs_8=float(max(abs(U[4]-U[2])))))
    if base!=4:
        for x in comparisons:
            for k in list(x):
                if k.endswith("_4_and_8"):
                    x[k[:-7]+"8_and_16"]=x.pop(k)
                elif k.endswith("_32_vs_8"):
                    x[k[:-7]+"32_vs_16"]=x.pop(k)
    return dict(T_GeV=T,layout=[f"q16/{base}",f"q32/{base}",f"q16/{2*base}",f"q32/{2*base}","q16/32"],cdf_momentum_over_T_thresholds=PROBES.tolist(),
        retained= data,comparisons=comparisons,input_sha256=[sha(p) for p in checkpoints],
        scope="Positive entropy-shape moments on actual matched early FRW warm restarts. All shape statistics are for retained particles only, not a full momentum continuum or a prediction for unavailable old exterior distributions.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("checkpoints",type=Path,nargs=5)
    p.add_argument("--base-blocks",type=int,default=4)
    a=p.parse_args()
    print(json.dumps(run(a.checkpoints,a.base_blocks),indent=2,allow_nan=False))
