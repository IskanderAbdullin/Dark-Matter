"""Vary how close to the start of reheating the dilute signed closure is initialized."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from archived_collisionless_streaming_reference import Background
from empty_lepton_pair_spectrum_reference import spectrum
from sm_dark_entropy_exchange_reference import MPL,StandardModel
from three_species_decay_onsager_reference import m,mp,M
from warm_prompt_reheating_expanding_slow_reference import visible,source_table
from warm_prompt_reheating_signed_slow_closure_reference import closed,quantum_state,pair_mb


def run(eos,phase_stop=.28):
    assert .075 < phase_stop <= .66
    bg=Background(eos)
    sm=StandardModel(eos)
    src=source_table(np.geomspace(.039,.101,67),1)
    masses=np.array((m,mp,M))
    origin=np.log(bg.scale(.1))
    stop_x=np.log(bg.scale(.04))-origin
    arms=[]
    for t0 in (.09998,.0999,.0995,.099,.097):
        seed=spectrum(sm,.1,t0,order=(12,35,24),nout=15)
        ee=np.asarray(seed["redshifted_means"])
        N0=np.array([sum(seed["stable_yields"]),seed["stable_yields"][1]])
        U0=ee[:2]@seed["stable_yields"]+ee[2]*N0[1]
        x0=np.log(bg.scale(t0))-origin
        def calc(x,Y,ret=False):
            T=float(np.exp(bg.v_to_u(origin+x)))
            rv,sv=visible(bg,T)
            Td,_,zm=closed(Y[:2]*sv,Y[2]*sv,T)
            Td,eta,N,E,P=quantum_state(Y[:2]*sv,Y[2]*sv,Td,zm,28 if phase_stop<.3 else 80,True)
            peaks=np.exp(eta-masses/Td)
            assert max(peaks)<max(.40,phase_stop+.06)
            J=np.exp(src(np.log(T)))
            H=np.sqrt((rv+E)/3)/MPL
            a=pair_mb(Td,42)*np.exp(2*np.array([eta[0],eta[1],eta[0]]))
            Delta=2*np.array([eta[2]-eta[0],eta[2]-eta[1],eta[1]-eta[0]])
            rates=-a*np.expm1(Delta)
            DY=2*np.array([rates[1]-rates[2],rates[0]+rates[2]])
            rhs=np.r_[(DY+J[:2])/(H*sv),J[2]/(H*sv)-3*P/sv]
            if ret:
                return T,Td,peaks,N/sv,DY/(H*sv),J[:2]/(H*sv)
            return rhs
        def stop(x,Y):
            return float(max(calc(x,Y,True)[2])-phase_stop)
        stop.terminal=True
        stop.direction=1
        def mark(x,Y):
            return float(max(calc(x,Y,True)[2])-.075)
        mark.direction=1
        MARK=[mark,stop]
        sol=solve_ivp(calc,(x0,stop_x),np.r_[N0,U0],method="DOP853",rtol=2e-8,atol=1e-16,max_step=.003,events=MARK,dense_output=True)
        assert sol.success,sol.message
        assert len(sol.t_events[0])==1
        if phase_stop==.28:
            assert len(sol.t_events[1])==1
        bounds=[float(a[0]) if len(a) else None for a in sol.t_events]
        checks=[]
        xs=np.unique(np.r_[np.linspace(x0,sol.t[-1],13 if phase_stop<.3 else 65),sol.t_events[0][0],sol.t[-1]])
        for j,Y in zip(xs,sol.sol(xs).T,strict=True):
            T,td,pk,pop,direct,dose=calc(j,Y,True)
            checks.append(dict(x_from_reheat=j,T_visible=T,T_dark=td,peak=pk.tolist(),physical_yields=pop.tolist(),
                               invariant_yields=Y[:2].tolist(),dark_and_visible_count_sources=[direct.tolist(),dose.tolist()]))
        arms.append(dict(anchor_t0=t0,anchor_passive_yields=N0.tolist(),
                        **{("x_075_and_028" if phase_stop==.28 else "x_075_and_selected"):bounds},rows=checks))
    return dict(scope="Controlled initialization-sensitivity test only. At prescribed early starts use the direct "
                      "FD lepton prompt spectra moments and impose a common dark quantum T and fast chemical affinity. "
                      "Then use the selected signed MB collision rates, FD lepton-only sources, adiabatic visible bath "
                      f"and dark pressure work. Stops at phase-factor {phase_stop}, if reached, or at 40 MeV. Not freeze-out or a relic. "
                      "These solved paths do not test the actual initial nonthermal relaxation or the missing hadronic, "
                      "elastic and in-medium physics.",arms=arms)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("--peak",type=float,default=.28)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.peak),indent=2,allow_nan=False))
