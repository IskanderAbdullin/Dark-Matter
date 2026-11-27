"""Re-evaluated below-pole elastic heat in the selected expanding cool closure."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator

from archived_collisionless_streaming_reference import Background
from sm_dark_entropy_exchange_reference import MPL
from warm_prompt_reheating_expanding_slow_reference import visible
from warm_prompt_reheating_quantum_heat_gap_reference import thermal
from warm_prompt_reheating_two_temperature_contact_probe_reference import fold
from warm_prompt_reheating_two_temperature_late_reference import fit,mb_both


def run(eos,prior):
    bg=Background(eos)
    b=next(a for a in json.loads(prior.read_text())["controls"][-1]["heat_arms"] if a["prompt_heat"])
    y=np.array([r["yields"]+r["energy_per_visible_entropy"] for r in b["rows"]])
    xb=np.array([r["x"] for r in b["rows"]])
    a0=np.log(bg.scale(b["rows"][0]["visible_T"]))
    coeffT=np.geomspace(.00014,.004,193)
    rate=PchipInterpolator(np.log(coeffT),np.log([mb_both(T,52) for T in coeffT]),axis=0,extrapolate=False)
    def decompose(x,u):
        Tv=np.exp(bg.v_to_u(a0+x))
        rv,s=visible(bg,Tv)
        l,p=[fit(u[j]*s,u[j+2]*s,Tv,j) for j in range(2)]
        H=np.sqrt((rv+s*u[2:].sum())/3)/MPL
        return Tv,s,l,p,H
    arms=[]
    f=PchipInterpolator(xb,y,axis=0)
    prev=None
    for count,qr in ((7,(7,10,10)),(11,(9,13,15)),(15,(9,13,15))):
        t=np.linspace(xb[0],xb[-1],count)
        cal=[]
        for x,u in zip(t,f(t),strict=True):
            tv,s,l,p,H=decompose(x,u)
            Tl,Tp=l[0],p[0]
            # At the common-temperature first sample infer the finite
            # conductivity from an explicitly tiny thermal tilt.
            pd=Tp if abs(Tp/Tl-1)>2e-6 else Tl*(1+2e-5)
            qq,ev,tail=fold(Tl,pd,[l[2],p[2]],qr,exact=True)
            z=-qq/((pd-Tl)*H*s)
            assert z>0 and ev>0.
            cal.append(dict(x=x,Tv=tv,Tl=Tl,Tp=Tp,heat_conductance_over_Hs=z,
                            fraction_of_scattering_proxy_above_cut=tail/(tail+ev)))
        kc=PchipInterpolator(t,np.log([r["heat_conductance_over_Hs"] for r in cal]))
        def func(x,u,ret=False):
            Tv,s,l,p,H=decompose(x,u)
            Tl,Tp=l[0],p[0]
            G=np.exp(rate(np.log(Tl))+2*l[1])-np.exp(rate(np.log(Tp))+2*p[1])
            assert np.all(np.isfinite(G))
            th=thermal(Tl,np.array([l[1],p[1]]),88,54,Tp/Tl-1)/(H*s)
            off=-np.exp(kc(x))*(Tp-Tl)
            g=G/(H*s)
            works=3*np.array([l[3],p[3]])/s
            d=np.r_[2*np.array([-g[0],g[0]]),np.array([-g[1]-th-off,g[1]+th+off])-works]
            return (Tv,Tl,Tp,d,g,th,off,works) if ret else d
        sol=solve_ivp(func,(xb[0],xb[-1]),y[0],method="DOP853",rtol=2e-8,atol=3e-16,max_step=.02,dense_output=True)
        assert sol.success and np.max(abs(sol.y[:2].sum(axis=0)/y[0,:2].sum()-1))<3e-10
        f=lambda x,dense=sol.sol:dense(x).T
        wdx=np.linspace(0.,xb[-1],21)
        tra=[]
        for x,u in zip(wdx,f(wdx),strict=True):
            Tv,Tl,Tp,d,g,th,off,works=func(x,u,True)
            tra.append(dict(x=x,visible_T=Tv,dark_T_light_scalar=[Tl,Tp],yields=u[:2].tolist(),
                            energy_per_visible_entropy=u[2:].tolist(),selected_pair_net_rate=g.tolist(),
                            prompt_heat=th,off_shell_scalar_heat=off,pressure_work=works.tolist()))
        v,w=np.polynomial.legendre.leggauss(4)
        grid=(sol.t[:-1,None]+sol.t[1:,None])/2+np.diff(sol.t)[:,None]*v/2
        w=np.diff(sol.t)[:,None]*w/2
        work=w.ravel()@np.array([sum(func(x,u,True)[-1]) for x,u in zip(grid.ravel(),f(grid.ravel()),strict=True)])
        err=abs(sol.y[2:,-1].sum()-y[0,2:].sum()+work)/y[0,2:].sum()
        assert err<2e-9
        if prev is not None:
            drift=max(abs(sol.y[:2,-1]/prev-1))
        else:drift=None
        prev=sol.y[:2,-1]
        arms.append(dict(path_samples=count,angle_and_radial=qr,calibration=cal,rows=tra,
                         relative_terminal_count_drift=drift,energy_work_defect=err))
    return dict(scope="Two-temperature thermal closure with iterated below-pole elastic heat and FD/BE occupation, "
                      "retaining only s<.003 GeV² for the added exact amplitude. Positive heat conductance is "
                      "recalculated on successive trajectories, not a Boltzmann evolution. Pair changes still use MB "
                      "moments and visible injection and heavier channels are omitted. No final relic claimed.",
                original_isolated_endpoint_yields=y[-1,:2].tolist(),arms=arms)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path);p.add_argument("prior",type=Path)
    a=p.parse_args();print(json.dumps(run(a.eos,a.prior),indent=2,allow_nan=False))
