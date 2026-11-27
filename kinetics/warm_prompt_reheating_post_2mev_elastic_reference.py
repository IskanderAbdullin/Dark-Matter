"""Bounded two-temperature continuation and separately tracked collision clocks."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

from archived_collisionless_streaming_reference import Background
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import m,mp
from warm_prompt_reheating_expanding_slow_reference import visible
from warm_prompt_reheating_two_temperature_contact_probe_reference import fold,kgeom
from warm_prompt_reheating_two_temperature_late_reference import fit,series,mb_both


def run(eos,prior):
    bg=Background(eos)
    x=json.loads(prior.read_text())["arms"][-1]["rows"][-1]
    t=x["visible_T"]
    a0=np.log(bg.scale(t))
    start=np.r_[x["yields"],x["energy_per_visible_entropy"]]
    y0=start[:2]
    end=np.log(bg.scale(.00008))-a0
    def get(x,u):
        tv=np.exp(bg.v_to_u(a0+x))
        rv,s=visible(bg,tv)
        l,p=[fit(u[i]*s,u[i+2]*s,tv,i) for i in range(2)]
        H=np.sqrt((rv+s*u[2:].sum())/3)/MPL
        return tv,s,l,p,H
    def guess(z):
        v=[]
        for j,T0 in enumerate(x["dark_T_light_scalar"]):
            T=T0*np.exp(-(2-j)*z)
            tv=np.exp(bg.v_to_u(a0+z))
            _,s=visible(bg,tv)
            sup=.899
            pk=brentq(lambda v:series(T,v,j)[0]/(y0[j]*s)-1,1e-10,sup,xtol=5e-15)
            n,K,_=series(T,pk,j)
            v.append((K+n*[m,mp][j])/s)
        return np.r_[y0,v]
    init=get(0.,start)
    def first(x):
        return np.array([start if q==0. else guess(q) for q in x])
    grid=np.geomspace(1e-7,.002,178)
    coeff=PchipInterpolator(np.log(grid),np.log([mb_both(T,48,True) for T in grid]),axis=0)
    out=[]
    prev=first
    # Preserve the same endpoint and compare these independent
    # per-trajectory rate calibrations rather than assuming fixed K.
    for n,order in ((8,(6,8,9)),(13,(8,10,13))):
        xs=np.linspace(0.,end,n)
        cal=[]
        for z,u in zip(xs,prev(xs),strict=True):
            Tv,s,l,p,H=get(z,u)
            TL,TS=l[0],p[0]
            ts=TS if abs(TS/TL-1)>1e-7 else TL*(1+3e-6)
            Q,N,O,K=fold(TL,ts,[l[2],p[2]],order,exact=True,transport=True)
            A=-Q/((ts-TL)*H*s)
            assert A>0.
            d=[]
            for j,(st,T,ma,g) in enumerate(((l,TL,m,4),(p,ts,mp,1))):
                k,E,w=kgeom(T,ma,12)
                f=st[2]*np.exp(-(E-ma)/T)
                w*=g*f/(1+f if j==0 else 1-f)
                d.append(w@(k*k))
            cal.append(dict(x=z,T_visible=Tv,Tl=TL,Ts=TS,heat_factor=A,
                            kick_over_H_second_moment=(K/(H*np.array(d))).tolist()))
        k=PchipInterpolator(xs,np.log([a["heat_factor"] for a in cal]),extrapolate=False)
        def calcr(x,u,report=False):
            TV,ss,l,p,H=get(x,u)
            TL,TS=l[0],p[0]
            vis=np.exp(coeff(np.log(TL))+2*np.log(l[2]))
            inv=np.exp(coeff(np.log(TS))+2*p[1]-2*m/TS)
            net=(vis-inv)/(H*ss)
            off=-np.exp(k(x))*(TS-TL)
            R=np.r_[2*np.array([-net[0],net[0]]),np.array([-net[1]-off,net[1]+off])-3*np.array([l[3],p[3]])/ss]
            return (TV,TL,TS,l[2],p[2],R,net,off) if report else R
        sol=solve_ivp(calcr,(0.,end),start,method="DOP853",rtol=5e-8,atol=2e-16,max_step=.045,dense_output=True)
        assert sol.success and max(abs(sol.y[:2].sum(axis=0)/sum(y0)-1))<1e-9
        prev=lambda v,dense=sol.sol:dense(v).T
        trail=[]
        for z,u in zip(np.linspace(0.,end,25),prev(np.linspace(0.,end,25)),strict=True):
            TV,TL,TS,fl,fs,r,G,Q=calcr(z,u,True)
            trail.append(dict(x=z,T_visible=TV,light_scalar_T=[TL,TS],phase_factors=[fl,fs],
                              yields=u[:2].tolist(),energy_per_visible_entropy=u[2:].tolist(),
                              net_per_H=G.tolist(),off_shell_heat=Q))
        out.append(dict(nodes=n,calibrations=cal,trajectory=trail))
    return dict(scope="Selected MB chemical moments with FD/BE thermal state and calibrated below-pole elastic heat. "
                      "No external bath, prompt exchange, heavier dark states, radiation feedback, or proof of internal "
                      "thermal distributions after momentum kicks slow. Small-temperature reaction coefficients factored "
                      "to avoid Boltzmann underflow. Ends at visible 80 keV for comparison only.",arms=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path);p.add_argument("handoff",type=Path)
    a=p.parse_args();print(json.dumps(run(a.eos,a.handoff),indent=2,allow_nan=False))
