"""A two-temperature controlled handoff after the fast heat clock slows."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import root
from scipy.special import expit,kve,roots_legendre

from archived_collisionless_streaming_reference import Background
from archived_scalar_pair_reaction_reference import amplitude
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import m,mp,M
from warm_prompt_reheating_expanding_slow_reference import visible,source_table
from warm_prompt_reheating_quantum_heat_gap_reference import thermal
from warm_prompt_reheating_signed_slow_closure_reference import pair_mb


def mb_both(T,n,scaled=False):
    w,wp=roots_legendre(n)
    x,xp=roots_legendre(n+16)
    ranges=np.array([0.,.5,2.,6.,12.,23.,42.,74.])
    d=np.diff(ranges)[:,None]/2
    u=2*m+T*((ranges[1:]+ranges[:-1])[:,None]/2+d*w)
    s=u*u
    p=np.sqrt(s/4-m*m)
    q=np.sqrt(s/4-mp*mp)
    invariant=4*(amplitude(s[:,:,None],x)@xp)*q/(64*np.pi*s*p)
    factor=2*u*T*d*wp*(s-4*m*m)*u*invariant/(32*np.pi**4)
    K0,K1=(kve(i,u/T)*np.exp(-(u-2*m if scaled else u)/T) for i in (0,1))
    return np.array([np.sum(factor*T*K1),np.sum(factor*(2*T*T*K1+T*u*K0))])


def series(T,peak,j):
    mass,g,sign=(m,4,-1.) if j==0 else (mp,1,1.)
    z=np.arange(1,113,dtype=float)
    a=z*mass/T
    n0=g*mass**2*T/z*kve(2,a)/(2*np.pi**2)
    n=n0*sign**(z-1)*np.exp(z*np.log(peak))
    e=3*T/z+mass*(kve(1,a)/kve(2,a)-1.)
    return np.array([n.sum(),n@e,n@(T/z)])


def fit(n,U,Tg,j):
    mass=m if j==0 else mp
    K=U-mass*n
    assert n>0 and K>0
    Tg=K/n/(1.7 if j==0 else 2.8)
    nu=series(Tg,.2,j)[0]
    peak=np.clip(.2*n/nu,1e-25,.78)
    def b(z):
        T=np.exp(z[0])
        f=expit(z[1])*.9
        return T,f,series(T,f,j)
    def residual(z):
        _,_,s=b(z)
        return np.log(s[:2]/[n,K])
    sol=root(residual,[np.log(Tg),np.log(peak/(.9-peak))],tol=3e-10)
    assert max(abs(residual(sol.x)))<3e-9
    T,p,k=b(sol.x)
    return T,float(np.log(p)+mass/T),p,k[2]


def run(eos,history):
    bg=Background(eos)
    d=json.loads(history.read_text())["controls"][-1]["rows"][-1]
    assert d["T_visible"]>0.005 and d["T_visible"]<0.007
    tv,T=d["T_visible"],d["dark_T"]
    eta=np.log(d["phase_factor"][:2])+np.array([m,mp])/T
    _,sv=visible(bg,tv)
    n=np.asarray(d["physical_yields"][:2])*sv
    assert d["physical_yields"][2]/sum(d["physical_yields"])<1e-8
    phases=np.array(d["phase_factor"][:2])
    initial=np.array([n/sv,[(series(T,phases[j],j)[1]+n[j]*ma)/sv for j,ma in enumerate((m,mp))]]).ravel()
    x0=np.log(bg.scale(tv))
    end=np.log(bg.scale(.002))-x0
    Ts=np.geomspace(.00014,.0037,187)
    out=[]
    cross=[]
    for t in (.00032,.001,.0025):
        f=mb_both(t,52)
        ref=pair_mb(t,52)[2]
        eps=1e-3
        derivative=(mb_both(1/(1/t-eps),52)[0]-mb_both(1/(1/t+eps),52)[0])/(2*eps)
        assert abs(f[0]/ref-1)<3e-10 and f[1]/f[0]>2*m and abs(derivative/f[1]-1)<3e-7
        cross.append([t,float(abs(f[0]/ref-1)),float(f[1]/f[0]),float(abs(derivative/f[1]-1))])
    src=source_table(np.geomspace(.00198,.00603,29),1)
    for nquad in (28,52):
        rate=PchipInterpolator(np.log(Ts),np.log([mb_both(t,nquad) for t in Ts]),axis=0,extrapolate=False)
        # Compare the same fixed clock with and without the exact prompt
        # heat. Neither handoff retains a visible source or heavy pairs.
        modes=[]
        for keep_heat in (False,True):
            def flux(x,u,ret=False):
                t=float(np.exp(bg.v_to_u(x0+x)))
                v,s=visible(bg,t)
                Y=u[:2]*s
                E=u[2:]*s
                l,p=[fit(Y[j],E[j],t,j) for j in range(2)]
                TL,TS=l[0],p[0]
                H=np.sqrt((v+E.sum())/3)/MPL
                G=np.exp(rate(np.log(TL))+2*l[1])-np.exp(rate(np.log(TS))+2*p[1])
                assert np.all(np.isfinite(G))
                h=thermal(TL,np.array([l[1],p[1]]),88,54,TS/TL-1) if keep_heat else 0.
                d=np.r_[2*np.array([-G[0],G[0]])/(H*s),
                      np.array([-G[1]-h,G[1]+h])/(H*s)-3*np.array([l[3],p[3]])/s]
                return (t,l,p,d,h/(H*s),G/(H*s)) if ret else d
            sol=solve_ivp(flux,(0.,end),initial,method="DOP853",rtol=3e-8,atol=1e-16,max_step=.014,dense_output=True)
            assert sol.success
            assert np.max(abs(sol.y[:2].sum(axis=0)/initial[:2].sum()-1))<2e-10
            rows=[]
            for x,u in zip(np.linspace(0.,end,17),sol.sol(np.linspace(0.,end,17)).T,strict=True):
                t,l,p,dv,h,G=flux(x,u,True)
                rows.append(dict(x=x,visible_T=t,light_and_scalar_dark_T=[l[0],p[0]],yields=u[:2].tolist(),
                                 light_and_scalar_phase_peak=[l[2],p[2]],energy_per_visible_entropy=u[2:].tolist(),
                                 net_pair_number_and_energy_over_Hs=G.tolist(),prompt_heat_over_Hs=h))
            ga,gw=roots_legendre(4)
            tnodes=(sol.t[1:]+sol.t[:-1])[:,None]/2+np.diff(sol.t)[:,None]*ga/2
            tw=np.diff(sol.t)[:,None]*gw/2
            projected=[]
            for xx,uu in zip(tnodes.ravel(),sol.sol(tnodes.ravel()).T,strict=True):
                t,l,s,dv,h,G=flux(xx,uu,True)
                vis,entropy=visible(bg,t)
                H=np.sqrt((vis+uu[2:].sum()*entropy)/3)/MPL
                J=np.exp(src(np.log(t)))
                projected.append([3*(l[3]+s[3])/entropy,*(J[:2]/(H*entropy)),J[2]/(H*entropy)])
            du=tw.ravel()@np.array(projected)
            defect=abs((sol.y[2:,-1].sum()-initial[2:].sum())+du[0])/initial[2:].sum()
            assert defect<1e-9
            modes.append(dict(prompt_heat=keep_heat,rows=rows,
                              conserved_minus_pressure_work_rel_defect=defect,
                              unblocked_leptonic_injection_fold_not_fed_back=du[1:].tolist()))
        out.append(dict(invariant_order=nquad,heat_arms=modes))
    return dict(invariant_number_and_mean_energy_checks=cross,scope="Selected two-temperature thermal model, with exact quantum thermodynamics for each species, "
                      "MB pair event rates in each direction at their own temperatures and optional exact FD/BE "
                      "prompt heat. No external visible inflow, heavy pair kinetics, or nonthermal evolution. "
                      "Ends at 2 MeV by construction, not a final relic or solved freeze-out.",controls=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("history",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.history),indent=2,allow_nan=False))
