"""Carry the isolated normal-phase moment model below its 40 MeV handoff."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator

from archived_collisionless_streaming_reference import Background
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import m,mp,M
from warm_prompt_reheating_expanding_slow_reference import visible,source_table
from warm_prompt_reheating_signed_slow_closure_reference import closed,quantum_state,pair_mb
from warm_prompt_reheating_saturation_distance_reference import compare
from warm_prompt_reheating_quantum_heat_gap_reference import thermal,capacity

MASSES=np.array([m,mp,M])


def initial(eos,data):
    B=Background(eos)
    p=data["arms"][0]["rows"][-1]
    T0=p["T_visible"]
    T=p["T_dark"]
    N=np.asarray(p["physical_yields"])
    _,s=visible(B,T0)
    z=np.asarray(p["peak"])*np.exp(MASSES/T)
    # Recover the starting dark energy from the independent species
    # densities and the resolved FD/BE energy integrals.
    from warm_prompt_reheating_saturation_distance_reference import measure
    from scipy.special import expit
    E,W,_=measure(T,38)
    eta=np.log(z)
    f=np.array([expit(eta[0]-E[0]/T),1/np.expm1(E[1]/T-eta[1]),expit(eta[2]-E[2]/T)])
    assert np.max(abs((W*f).sum(axis=(1,2))/s/N-1))<1e-8
    return p,T0,T,z,float((W*f*E).sum()/s)


def run(eos,history,cold=False):
    B=Background(eos)
    p,T0,T0d,z0,U0=initial(eos,json.loads(history.read_text()))
    a0=np.log(B.scale(T0))
    y0=np.r_[p["invariant_yields"],U0]
    floor=.002 if cold else .008
    end=np.log(B.scale(floor))-a0
    src=source_table(np.geomspace(.98*floor,.0403,73 if cold else 59),1)
    dgrid=np.geomspace(.00035,.012,221)
    cores=[]
    for n in (26,52):
        coeff=np.array([pair_mb(z,n) for z in dgrid])
        assert np.min(coeff)>0.
        mtable=PchipInterpolator(np.log(dgrid),np.log(coeff),extrapolate=False,axis=0)
        def rhs(x,Y,ret=False):
            T=float(np.exp(B.v_to_u(x+a0)))
            rv,sv=visible(B,T)
            TM,_,zM=closed(Y[:2]*sv,Y[2]*sv,T)
            Td,eta,P,E,pressure=quantum_state(Y[:2]*sv,Y[2]*sv,TM,zM,112,True)
            peak=np.exp(eta-MASSES/Td)
            # Do not solve at or past a saturated normal-mode occupation.
            assert max(peak)<.88
            H=np.sqrt((rv+E)/3)/MPL
            J=np.exp(src(np.log(T)))
            K=np.exp(mtable(np.log(Td)))
            assert np.all(np.isfinite(K)) and np.all(K>0.)
            f=K*np.exp(2*np.array([eta[0],eta[1],eta[0]]))
            r=K*np.exp(2*np.array([eta[2],eta[2],eta[1]]))
            G=f-r
            dv=2*np.array([G[1]-G[2],G[0]+G[2]])
            rate=np.r_[(dv+J[:2])/(H*sv),J[2]/(H*sv)-3*pressure/sv]
            if ret:
                heat=None
                if cold:
                    Kha=thermal(Td,eta,68,45)
                    caph=np.array([capacity(Td,eta[0],m,4,False,20),capacity(Td,eta[1],mp,1,True,20)])
                    mode=Kha*np.sum(1/caph)/H
                    heat=(float(mode),float(mode/max(1.,max(abs(rate[:2]/Y[:2])))))
                return T,Td,P/sv,E/sv,sv,peak,rate,(f/(H*sv)).tolist(),(r/(H*sv)).tolist(),heat
            return rate
        def stop(x,Y): return float(max(rhs(x,Y,True)[5])-.80)
        stop.terminal=True
        stop.direction=1
        def heat_stop(x,Y):
            return rhs(x,Y,True)[-1][1]-30.
        heat_stop.terminal=True
        heat_stop.direction=-1
        events=(stop,heat_stop) if cold else (stop,)
        sol=solve_ivp(rhs,(0.,end),y0,method="DOP853",rtol=2e-8,atol=2e-16,max_step=.012,events=events,dense_output=True)
        assert sol.success
        out=[]
        for index,x in enumerate(np.linspace(0.,sol.t[-1],29)):
            y=sol.sol(x)
            T,Td,phys,E,sv,peak,delta,forward,reverse,heat=rhs(x,y,True)
            rec=dict(x=x,T_visible=T,dark_T=Td,physical_yields=phys.tolist(),invariants=y[:2].tolist(),
                      phase_factor=peak.tolist(),net_invariant_rates=delta[:2].tolist(),pair_forward=forward,pair_reverse=reverse)
            if cold:rec["heat_over_H_and_over_max_H_or_other_rates"]=list(heat)
            if index in (0,14,28):
                rec["thermal_capacity_test"]=compare(Td,y[0]*sv,y[1]*sv,E*sv)
            out.append(rec)
        cores.append(dict(order=n,stopped_at_phase_factor=bool(len(sol.t_events[0])),
                          stopped_at_heat_clock_cut=bool(len(sol.t_events[1])) if cold else False,rows=out))
    return dict(scope="Continuation of earlier MB-event-rates + FD/BE thermodynamic closure, beginning with the "
                      "actual modeled near-empty arm at 40 MeV. Pair kernels are interpolated with two resolutions. "
                      "The initial state inherits all earlier modeling omissions. Not an exact collision evolution "
                      "or final relic; cold mode also stops at a declared heat-to-other-rates ratio of 30. "
                      "A declared clock ratio is not a physical rate crossing at unity or a freeze-out.",
                controls=cores)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("history",type=Path)
    p.add_argument("--cold",action="store_true")
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.history,a.cold),indent=2,allow_nan=False))
