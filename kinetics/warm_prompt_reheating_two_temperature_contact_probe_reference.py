"""Off-shell low-invariant isotropic contact heat diagnostic on the cool two-T path."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.special import expit,roots_legendre

from archived_collisionless_streaming_reference import Background
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import m,mp,M,y
from warm_prompt_reheating_expanding_slow_reference import visible
from warm_prompt_reheating_two_temperature_late_reference import series


def exact_elastic_sigma_at_angle(S,z):
    """Spin-averaged 4π-weighted dσ/dΩ using the true scalar scattering cosine."""
    pc2=np.maximum((S-(m+mp)**2)*(S-(m-mp)**2)/(4*S),0.)
    t=-2*pc2*(1-z)
    u=2*(m*m+mp*mp)-S-t
    d1,d2=S-M*M,u-M*M
    a=M*(1/d1+1/d2)
    b=(S+m*m-mp*mp)/(2*d1)+(u+m*m-mp*mp)/(2*d2)
    Q2=S/d1**2+u/d2**2+2*(m*m-mp*mp)/(d1*d2)
    dot=m*m-t/2
    trace=2*a*a*(dot+m*m)+8*a*m*b+2*(2*b*b+(m*m-dot)*Q2)
    assert np.min(trace)>=0.
    return y**4*trace/(16*np.pi*S)


def kgeom(T,mass,n):
    t,w=roots_legendre(n)
    e=np.r_[0.,np.geomspace(.003,38.,17)]
    h=np.diff(e)/2
    k=T*((e[1:]+e[:-1])[:,None]/2+h[:,None]*t)
    p=np.sqrt(k*(2*mass+k))
    # p² dp = p E dE.
    dn=w*h[:,None]*T*p*(k+mass)/(2*np.pi**2)
    return p.ravel(),(k+mass).ravel(),dn.ravel()


def fold(T,U,pk,order,exact=False,smax=.003,azimuth=10):
    p,E,w=kgeom(T,m,order[0]);q,G,v=kgeom(U,mp,order[0])
    x,wx=roots_legendre(order[1]);z,wz=roots_legendre(order[2])
    f=lambda en:expit(np.log(pk[0])-(en-m)/T)
    g=lambda en:(lambda a:a/(1-a))(pk[1]*np.exp(-(en-mp)/U))
    S0=y**4*m*m/(np.pi*(M-m)**2*(m+mp)**2)
    count,heat,countout=0.,0.,0.
    for a in range(0,len(p),8):
        L=p[a:a+8,None,None]
        U0=E[a:a+8,None,None]
        R=q[None,:,None];U1=G[None,:,None]
        dot=U0*U1-L*R*x
        S=m*m+mp*mp+2*dot
        Q=U0+U1
        P=np.sqrt(np.maximum(Q*Q-S,0.))
        cms=(S+mp*mp-m*m)/(2*np.sqrt(S))
        pcm=np.sqrt(np.maximum(cms*cms-mp*mp,0.))
        eu=Q*cms/np.sqrt(S)
        span=P*pcm/np.sqrt(S)
        out=np.maximum(eu[:,:,:,None]+span[:,:,:,None]*z,mp)
        dep=out-U1[:,:,:,None]
        elast1=f(U0)[:,:,:,None]*g(U1)[:,:,:,None]*(1-f(Q[:,:,:,None]-out))*(1+g(out))
        elast2=f(Q[:,:,:,None]-out)*g(out)*(1-f(U0))[:,:,:,None]*(1+g(U1))[:,:,:,None]
        vel=np.sqrt(np.maximum(dot*dot-(m*mp)**2,0.))/(U0*U1)
        step=(4*w[a:a+8,None,None]*v[None,:,None]*wx/2)*vel
        inside=S<smax
        coef=S0
        if exact:
            ci=np.clip(np.divide(U1-eu,span,out=np.zeros_like(eu),where=span>0),-1.,1.)
            ci_z=ci[:,:,:,None]*z
            perp=np.sqrt(np.maximum(1-ci*ci,0.))[:,:,:,None]*np.sqrt(1-z*z)
            Ss=np.where(inside,S,.002)[:,:,:,None]
            coeff=np.zeros_like(ci_z)
            for k in np.cos((np.arange(azimuth)+.5)*np.pi/azimuth):
                coeff+=exact_elastic_sigma_at_angle(Ss,np.clip(ci_z+perp*k,-1.,1.))/azimuth
            coef=coeff
        w0=(step*inside)[:,:,:,None]*wz/2*coef
        heat+=.5*np.sum(w0*(elast1-elast2)*dep)
        count+=np.sum(w0*elast1)
        countout+=np.sum((step*~inside)[:,:,:,None]*wz/2*S0*elast1)
    return float(heat),float(count),float(countout)


def run(eos,prior,exact=False):
    b=Background(eos)
    d=json.loads(prior.read_text())
    r=next(x for x in d["controls"][-1]["heat_arms"] if x["prompt_heat"])["rows"]
    out=[]
    for i in ((4,8,16) if exact else (0,4,8,12,16)):
        v=r[i]
        tv=v["visible_T"]
        tl,ts=v["light_and_scalar_dark_T"]
        pk=v["light_and_scalar_phase_peak"]
        rr,s=visible(b,tv)
        H=np.sqrt((rr+s*sum(v["energy_per_visible_entropy"]))/3)/MPL
        levels=[]
        for n in (((8,10,12),(12,15,17)) if exact else ((10,12,13),(15,18,21))):
            Q,W,X=fold(tl,ts,pk,n,exact)
            assert Q*(tl-ts)>=-1e-45
            levels.append(dict(orders=n,heat_Q_over_Hs=Q/(H*s),
                               heat_over_scalar_energy_per_H=Q/(H*s*v["energy_per_visible_entropy"][1]),
                               gross_events_over_scalar_count_per_H=W/(H*s*v["yields"][1]),
                               full_contact_count_fraction_outside_s_cut=X/(W+X)))
        out.append(dict(T_visible=tv,T_light=tl,T_scalar=ts,controls=levels))
    tag=("Invariant differential s/u amplitude below the cut, with explicit boost and azimuth."
         if exact else "Isotropic threshold total cross section (contact control).")
    return dict(exact_unpolarized_amplitude_below_cut=exact,scope=tag+" Exact lab kinematics and FD/BE occupation on held "
                "two-temperature states. s<.003 GeV² excludes the on-shell pole. The out-of-cut contact proxy "
                "is not a bound on omitted interactions.",rows=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("prior",type=Path)
    p.add_argument("--exact",action="store_true")
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.prior,a.exact),indent=2,allow_nan=False))
