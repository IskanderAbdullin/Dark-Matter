"""Quantum moment matching with alternate local occupation-factor choices."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from archived_collisionless_streaming_reference import Background
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import m,mp,M
from warm_prompt_reheating_expanding_slow_reference import visible,source_table
from warm_prompt_reheating_signed_slow_closure_reference import closed,quantum_state,pair_mb


def one(bg,initial,scenario):
    t0=initial["TF"]
    a0=np.log(bg.scale(t0))
    x_end=np.log(bg.scale(.035))-a0
    y=np.array([sum(initial["stable_yields"]),initial["stable_yields"][1]])
    E=np.asarray(initial["redshifted_means"])
    e=E[:2]@initial["stable_yields"]+E[2]*y[1]
    src=source_table(np.geomspace(.0345,.0605,39),1)

    def rhs(x,Y,ret=False):
        t=float(np.exp(bg.v_to_u(a0+x)))
        rv,ss=visible(bg,t)
        tm,_,zm=closed(Y[:2]*ss,Y[2]*ss,t)
        T,eta,n,_,P=quantum_state(Y[:2]*ss,Y[2]*ss,tm,zm,16,True)
        z=np.exp(eta)
        pk=z*np.exp(-np.array([m,mp,M])/T)
        assert pk.max()<.15
        l,s,h=pk
        low=np.array([1/((1+l)*(1+h))**2,1/(1+h)**2,1/(1+l)**2])
        high=np.array([1.,1/(1-s)**2,1/(1-s)**2])
        omega=low+(high-low)*scenario
        mb=pair_mb(T,52)*np.array([z[0]**2,z[1]**2,z[0]**2])
        detail=-np.expm1(2*np.array([eta[2]-eta[0],eta[2]-eta[1],eta[1]-eta[0]]))
        W=mb*detail*omega
        H=np.sqrt((rv+Y[2]*ss)/3)/MPL
        J=np.exp(src(np.log(t)))
        rates=np.r_[(2*np.array([W[1]-W[2],W[0]+W[2]])+J[:2])/(H*ss),J[2]/(H*ss)-3*P/ss]
        return (rates,t,T,pk,n/ss,high,low) if ret else rates

    def stop(x,Y):
        return max(rhs(x,Y,True)[3])-.075
    stop.terminal=True
    stop.direction=1
    sol=solve_ivp(rhs,(0.,x_end),np.r_[y,e],events=stop,method="DOP853",rtol=2e-9,atol=1e-16,max_step=.002,dense_output=True)
    assert sol.success and len(sol.t_events[0])==1
    rows=[]
    for x,v in zip(np.linspace(0,sol.t[-1],11),sol.sol(np.linspace(0,sol.t[-1],11)).T,strict=True):
        r,t,T,pk,n,hi,lo=rhs(x,v,True)
        rows.append(dict(log_scale_growth=x,T_visible_GeV=t,T_dark_GeV=T,max_phase_occupation=pk.tolist(),
                         invariant_light_scalar_yields=v[:2].tolist(),actual_light_scalar_heavy_yields=n.tolist(),
                         signed_rate=r[:2].tolist(),occupation_factor_lower_upper=[lo.tolist(),hi.tolist()]))
    return dict(local_forward_factor_position=scenario,rows=rows)


def run(eos,seed):
    b=Background(eos)
    x=json.loads(seed.read_text())
    return dict(scope="Thermodynamic closure uses dilute-series FD/BE moments and pressure. Pair rates use "
                      "both directions through exact thermal affinity, with individual forward factors "
                      "chosen at fixed positions in locally valid occupation bounds. These three model "
                      "trajectories are rate-sensitivity probes, not upper/lower bounds on a solution, "
                      "and are not an exact kinetic or quantum collision integration.",
                paths=[one(b,x,q) for q in (0.,.5,1.)])


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path);p.add_argument("seed",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.seed),indent=2,allow_nan=False))
