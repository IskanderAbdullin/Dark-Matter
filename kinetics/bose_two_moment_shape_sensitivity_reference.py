"""Positive Bose diagnostic closures with identical N and kinetic E."""
import json
import numpy as np
from scipy.optimize import brentq
from scipy.special import expit, logsumexp, roots_legendre
from massive_fv_fixed_bath_reference import bounds, kinematics


def qmesh(lo,hi,m,order):
    z,w=roots_legendre(order)
    q=lo+(hi-lo)*(z+1)/2
    dw=(hi-lo)/2*w
    e=np.hypot(q,m)
    k=q*q/(e+m)
    kl=lo*lo/(np.hypot(lo,m)+m) if lo else 0.
    kh=hi*hi/(np.hypot(hi,m)+m)
    return q,dw,e,k,(k-kl)/(kh-kl),kl,kh


def infer(lo,hi,m,fill,mean_s,thermal,order):
    q,dw,e,_,s,kl,kh=qmesh(lo,hi,m,order)
    wq=q*q*dw
    norm=wq.sum()
    prior=-np.log(np.expm1(e)) if thermal else np.zeros_like(e)
    l=np.log(wq)+prior
    def moment(beta):
        u=l-beta*s
        return float(np.exp(u-logsumexp(u))@s)
    a,b=-1.,1.
    for _ in range(20):
        if moment(a)>=mean_s>=moment(b): break
        a,b=2*a,2*b
    else:
        raise ValueError("energy target outside quadrature-node support")
    beta=brentq(lambda r:moment(r)-mean_s,a,b,xtol=1e-12)
    A=np.log(fill*norm)-logsumexp(l-beta*s)
    def f(e,k):
        logp=-np.log(np.expm1(e)) if thermal else 0.
        return np.exp(A+logp-beta*(k-kl)/(kh-kl))
    return f,beta


def control(lo,hi,m,M,b,fill,target,order):
    ans=[]
    q,dw,e,k,ss,_,_=qmesh(lo,hi,m,2*order)
    K,D=kinematics(M,m,b)
    r_l,r_h=bounds(e,m,b,K,D)
    z,w=roots_legendre(2*order)
    r=r_l[:,None]+(r_h-r_l)[:,None]*(z+1)/2
    rwt=(r_h-r_l)[:,None]*w/2
    G=np.sum(rwt*expit(-e[:,None]-r)*expit(r),axis=1)
    qe=dw*q/e*G
    pref=4/(2*np.pi**2)*M/(2*D/M**2)
    for thermal in (False,True):
        f,slope=infer(lo,hi,m,fill,target,thermal,order)
        values=f(e,k)
        phase=dw*q*q
        nmean=float(phase@values/phase.sum())
        emean=float((phase*ss)@values/(phase@values))
        assert np.all(values>0) and abs(nmean/fill-1)<5e-10
        assert abs(emean-target)<5e-10
        rate=float(pref*(qe@(1-np.expm1(e)*values)))
        ans.append(dict(thermal_reference=thermal,
                        cell_mean=nmean,kinetic_fraction=emean,slope=slope,number_rate=rate))
    diff=abs(ans[0]["number_rate"]-ans[1]["number_rate"])/max(abs(t["number_rate"]) for t in ans)
    return dict(qlo=lo,qhi=hi,fill=fill,mean_s=target,order=order,closures=ans,
                relative_rate_difference_over_larger=diff)


def run():
    M,m,b=.005,.0004,.0009
    pars=[(0.,.008,250.,.35),(.0001,.04,80.,.55),(.005,.15,50.,.6)]
    a=[control(lo,hi,m,M,b,n,e,64) for lo,hi,n,e in pars]
    c=[control(lo,hi,m,M,b,n,e,96) for lo,hi,n,e in pars]
    drift=max(abs(v["number_rate"]/w["number_rate"]-1) for x,y in zip(a,c)
              for v,w in zip(x["closures"],y["closures"]))
    assert drift<2e-7
    assert max(z["relative_rate_difference_over_larger"] for z in c)>.001
    return dict(rows=c,max_relative_quadrature_drift=drift,
                note="Both priors are diagnostics only; no physical freeze-out claim.")

if __name__ == "__main__":
    print(json.dumps(run(),indent=2,allow_nan=False))
