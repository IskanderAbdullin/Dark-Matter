"""Path-tabulated selected quantum rate corrections for a scoped thermal closure."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.special import expit

from archived_collisionless_streaming_reference import Background
from archived_scalar_pair_reaction_reference import gross,points
from sm_dark_entropy_exchange_reference import MPL,StandardModel
from empty_lepton_pair_spectrum_reference import spectrum
from three_species_decay_onsager_reference import m,mp,M
from warm_prompt_reheating_expanding_slow_reference import visible,source_table
from warm_prompt_reheating_heavy_gain_probe_reference import rates
from warm_prompt_reheating_signed_slow_closure_reference import closed,quantum_state,pair_mb

MASSES=np.array([m,mp,M])
SCALAR_LIMIT=.28


def factors(T,eta,index):
    peaks=np.exp(eta-MASSES/T)
    assert peaks.max()<.78
    l,s,h=peaks
    F=lambda E:1./np.expm1(E/T-eta[1])
    G=lambda E:expit(eta[0]-E/T)
    n,rad,angle,azi=((54,3,17,12),(74,4,25,16))[index]
    edge=np.r_[0.,np.geomspace(.00035,max(65*T,2.),n)]
    pp,_,np0=points(edge,m,G,rad)
    ss,_,ns0=points(edge,mp,F,rad)
    ar=np.r_[rates(pp,4*np0,m,angle).sum(),rates(ss,ns0,mp,angle).sum()]
    _,q,_=gross(edge,m,mp,G,F,False,rad,angle,azi,9,False,exchange_mass=M)
    # Independent reverse light/scalar check. At common T the logarithmic
    # detailed-balance identity does not depend on quadrature or phase.
    _,v,_=gross(edge,mp,m,F,G,True,rad,angle,azi,9,False,exchange_mass=M)
    err=abs(v[1]/q[1]/np.exp(2*(eta[1]-eta[0]))-1)
    assert err<.003
    z=np.exp(eta)
    bare=pair_mb(T,60)*np.array((z[0]**2,z[1]**2,z[0]**2))
    w=np.r_[ar*((1+(1+h)**-2)/2),q[1]]/bare
    lo=np.array([(1+l)**-2*(1+h)**-2,(1+h)**-2,(1+l)**-2])
    hi=np.array([1.,(1-s)**-2,(1-s)**-2])
    # Finite excited-state shell quadrature can slightly overshoot an
    # occupation upper. Record and test against both a bound and tolerance.
    defect=float(max(0.,max(lo-w),max(w-hi)))
    assert defect<.012
    return dict(weights=w.tolist(),lo=lo.tolist(),hi=hi.tolist(),detailed_balance_relative_error=float(err),
                approximate_heavy_Pauli_factor_total_width=1-(1+h)**-2,max_finite_quad_violation=defect)


def build(bg,seed,wide=False):
    t0=seed["TF"]
    a0=np.log(bg.scale(t0))
    y=np.array([sum(seed["stable_yields"]),seed["stable_yields"][1]])
    ms=seed["redshifted_means"]
    energy=np.array(ms[:2])@seed["stable_yields"]+ms[2]*y[1]
    src=source_table(np.geomspace(.0345,t0+.0005,32 if t0<.07 else 65),1)
    def unpack(x,Y):
        t=float(np.exp(bg.v_to_u(a0+x)))
        rv,sv=visible(bg,t)
        tm,_,zm=closed(Y[:2]*sv,Y[2]*sv,t)
        T,eta,N,E,P=quantum_state(Y[:2]*sv,Y[2]*sv,tm,zm,80 if wide else 28,True)
        z=np.exp(eta)
        pk=z*np.exp(-MASSES/T)
        assert pk.max()<(.78 if wide else .45)
        return t,rv,sv,T,eta,N,E,P,pk
    return t0,a0,np.r_[y,energy],src,unpack


def run(eos,seed,early=False,wide=False):
    assert not wide or early
    bg=Background(eos)
    if early:
        prepared=spectrum(StandardModel(eos),.1,.09998,order=(12,35,24),nout=15)
    else:
        prepared=json.loads(seed.read_text())
    t0,a0,y0,src,unpack=build(bg,prepared,wide)
    end=np.log(bg.scale(.04 if wide else .035))-a0
    limit=.65 if wide else SCALAR_LIMIT
    def integrate(table,order):
        def fun(x,Y,ret=False):
            t,rv,sv,T,eta,N,E,P,pk=unpack(x,Y)
            z=np.exp(eta)
            H=np.sqrt((rv+E)/3)/MPL
            J=np.exp(src(np.log(t)))
            mb=pair_mb(T,order)*np.exp(2*np.array([eta[0],eta[1],eta[0]]))
            net=mb*(-np.expm1(2*np.array([eta[2]-eta[0],eta[2]-eta[1],eta[1]-eta[0]])))
            if table is not None:
                x0,x1=table.x[[0,-1]]
                net*=table(np.clip(x,x0,x1))
            d=np.r_[(2*np.array([net[1]-net[2],net[0]+net[2]])+J[:2])/(H*sv),J[2]/(H*sv)-3*P/sv]
            return (t,T,N/sv,pk,d) if ret else d
        def active(x,Y):
            return max(unpack(x,Y)[-1])-limit
        active.terminal=True
        active.direction=1
        s=solve_ivp(fun,(0.,end),y0,rtol=3e-9,atol=2e-16,method="DOP853",
                    max_step=.002,dense_output=True,events=active)
        assert s.success and (wide or len(s.t_events[0])==1)
        return s,fun
    source,ff=integrate(None,52)
    history=[]
    # Evaluate a new rate table on the recalculated trajectory and repeat
    # with separate numbers of knots and event-angle quadrature.
    for grid,index in (((11,0),(16,1)) if wide else (((9,0),(13,1)) if early else ((9,0),(13,1),(17,1)))):
        nodes=np.linspace(0.,source.t[-1],grid)
        rows=[]
        for x in nodes:
            Y=source.sol(x)
            t,rv,sv,T,eta,N,E,P,pk=unpack(x,Y)
            rows.append(dict(log_scale_growth=x,T_dark=T,max_classical_phase_factor=pk.tolist(),**factors(T,eta,index)))
        table=PchipInterpolator(nodes,[q["weights"] for q in rows],axis=0)
        source,ff=integrate(table,52)
        trace=[]
        for x in np.linspace(0.,source.t[-1],65 if wide else 15):
            Y=source.sol(x)
            t,T,N,pk,d=ff(x,Y,True)
            trace.append(dict(log_scale_growth=x,visible_T=t,dark_T=T,invariants=Y[:2].tolist(),actual_counts=N.tolist(),
                              phase_peak=pk.tolist(),d_invariants_over_Hs=d[:2].tolist()))
        history.append(dict(grid=grid,angle_high=bool(index),table=rows,samples=trace))
    return dict(phase_limit=limit,wide=wide,scope="Thermal FD/BE moment closure, exact thermal affinities, selected vacuum rates with tabulated "
                      "quantum multipliers iterated on their own path. LS uses direct in/out occupation integration. "
                      "Heavy production uses direct incoming distributions and midpoint of its final Pauli envelope. "
                      "No dark/visible elastic heating, fully coupled visible cooling, hadronic injection, or "
                      "validation of an actual thermal kinetic distribution. The chosen phase-factor stop or minimum temperature is a "
                      "test boundary, not a physical transition, freeze-out, or a final relic yield.",history=history)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("seed",type=Path,nargs="?")
    p.add_argument("--early",action="store_true")
    p.add_argument("--wide",action="store_true")
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.seed,a.early,a.wide),indent=2,allow_nan=False))
