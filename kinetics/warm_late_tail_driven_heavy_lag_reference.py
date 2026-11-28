"""Lag of the heavy parent driven by the saved short late tail trajectory."""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from joint_visible_heavy_adiabatic_reference import angle_panels, snapshot_functions
from late_radiation_background_reference import MPL, TabulatedRadiation
from three_species_decay_onsager_reference import Es, M, P, m, mp, width
from warm_coupled_serial_restart_reference import atomic, sha
from warm_light_occupation_stage_audit_reference import digest, read
from warm_prompt_reentry_shape_audit_reference import get_mom


def weighted_update(d0,h0,h1,g0,g1,L):
    """Endpoint lag and mean lag times linear g for one prescribed interval."""
    regular=L>=.025
    small=~regular
    inv=np.divide(1.,L,out=np.zeros_like(L),where=regular)
    e=np.exp(-L)
    k0=-np.expm1(-L)*inv
    k1=(k0-e)*inv
    s0=(1.-k0)*inv
    s1=(.5-k1)*inv
    x=L[small]
    if np.any(small):
        # Taylor through the eighth-order term, finite also at zero.
        from math import factorial
        series=lambda m,k: np.polynomial.polynomial.polyval(
            x,[(-1.)**i/(factorial(i+k)*(i+m+k)) for i in range(9)])
        k0[small]=series(1,0)
        k1[small]=series(2,0)
        s0[small]=series(1,1)
        s1[small]=series(2,1)
    delta=h1-h0
    d1=d0*e-delta*k0
    mean=d0*(g0*k0+(g1-g0)*k1)-delta*(g0*s0+(g1-g0)*s1)
    assert np.all((L>0.)|((g0==0.)&(g1==0.)&(d0==0.)&(delta==0.)))
    return d1,mean


def check():
    from scipy.integrate import solve_ivp
    for L in (1e-10,.021,.03,3.2,350.):
        d,h0,h1,g0,g1=.003,.025,.041,.6,1.3
        z=solve_ivp(lambda t,y:[-L*y[0]-(h1-h0),(g0+(g1-g0)*t)*y[0]],
                    (0.,1.),[d,0.],method="Radau",atol=2e-15,rtol=2e-13)
        E,W=weighted_update(*map(np.atleast_1d,(d,h0,h1,g0,g1,L)))
        np.testing.assert_allclose([E[0],W[0]],z.y[:,-1],atol=2e-14,rtol=1e-11)


def geometry(ep,el,cut,order):
    A=np.unique(np.r_[ep,cut])
    B=np.unique(np.r_[el,cut])
    def gate(wall,mass,rest):
        p=np.sqrt((wall-mass)*(wall+mass))
        u=wall*rest+p*P
        return np.r_[M*(rest*rest+p*p)/u,M*u/mass**2]
    lo,hi=max(cut,M),ep[-1]+el[-1]
    walls=np.unique(np.clip(np.r_[lo,hi,gate(A,mp,Es),gate(B,m,M-Es),
                                   (A[:,None]+B).ravel()],lo,hi))
    z,w=roots_legendre(order)
    d=np.diff(walls)/2
    v=((walls[:-1,None]+walls[1:,None])/2+d[:,None]*z).ravel()
    dv=(d[:,None]*w).ravel()
    measure=dv*2/np.pi**2*np.sqrt((v-M)*(v+M))*width*M
    return v,measure,A,B


def rates(states,v,A,B,ep,el,cut,order):
    out=[np.empty((4,len(v))) for _ in states]
    for i in range(0,len(v),120):
        E=v[i:i+120]
        phi,W=angle_panels(E,A,B,angle=order)
        chi=E[:,None,None]-phi
        W*=((phi>=ep[0])&(phi<=ep[-1])&(chi>=el[0])&(chi<=el[-1]))
        excess=np.maximum(phi-cut,0.)+np.maximum(chi-cut,0.)
        for (f,g),dest in zip(states,out,strict=True):
            b,c=f(phi),g(chi)
            ab=W*b*c
            loss=W*(1+b)*(1-c)
            full=ab+loss
            n=full.sum(axis=(1,2))
            h=np.divide(ab.sum(axis=(1,2)),n,out=np.zeros_like(n),where=n>0.)
            gw=(full*excess).sum(axis=(1,2))
            stat=((h[:,None,None]*loss-(1-h[:,None,None])*ab)*excess).sum(axis=(1,2))
            dest[:,i:i+len(E)]=[h,n,gw,stat]
    return out


def run(eos,paths,order,pieces,initial):
    check()
    raw=[read(p) for p in paths]
    assert all(x["hash"]==digest(x["snapshot"]) and x["plan"]==raw[0]["plan"] for x in raw)
    assert raw[0]["plan"]["eos_sha256"]==sha(eos)
    S=raw[0]["snapshot"]
    t,grid=S["T_visible_GeV"],np.array(S["lab_momentum_edges_GeV"])
    C=raw[0]["plan"]["normalized_cut_GeV"]
    keys=("entropy_variable_scalar","entropy_variable_fermion")
    phase=[np.asarray([x["snapshot"][k] for k in keys]) for x in raw]
    assert all(x["snapshot"]["T_visible_GeV"]==t and np.array_equal(x["snapshot"]["lab_momentum_edges_GeV"],grid) for x in raw)
    rho=get_mom(grid,phase[0])[2].sum()
    H=np.sqrt((TabulatedRadiation(eos,t).state(1.)[1]+rho)/3)/MPL
    _,ep,_,el=snapshot_functions(S)
    v,measure,A,B=geometry(ep,el,C,order)
    states=[]
    for a,b in zip(phase,phase[1:]):
        for j in range(pieces):
            z=j/pieces
            t0,t1=(1.-z)*a+z*b
            f,_,g,_=snapshot_functions(dict(S,**{keys[0]:t0,keys[1]:t1}))
            states.append((f,g))
    f,_,g,_=snapshot_functions(raw[-1]["snapshot"])
    states.append((f,g))
    arr=rates(states,v,A,B,ep,el,C,order)
    a=arr[0]
    lag=np.zeros_like(v) if initial=="equilibrium" else -a[0].copy()
    moments=[]
    Emeasure=measure*v*v/(width*M*rho)
    rate_factor=measure/(rho*H)
    current=width*M/(v*H)
    cumulative=0.
    for k,(i,j) in enumerate(zip(raw,raw[1:])):
        step=(j["time"]-i["time"])/pieces
        assert step>0.
        sums=np.zeros(2)
        for s in range(pieces):
            b=arr[k*pieces+s+1]
            L=current*np.minimum(a[1],b[1])*step
            lag,weighted=weighted_update(lag,a[0],b[0],a[2],b[2],L)
            sums+=step*np.array([rate_factor@weighted,rate_factor@((a[3]+b[3])/2)])
            a=b
        cumulative+=sums[0]
        moments.append(dict(time=j["time"],minimum_response_per_substep=float(np.min(L[L>0.])),
              signed_W_correction_at_endpoint=float(rate_factor@(lag*b[2])),
              signed_prompt_at_endpoint=float(rate_factor@b[3]),
              W_correction_integrated_over_rho=float(sums[0]),
              W_stationary_integrated_over_rho=float(sums[1]),
              cumulative_W_correction_over_rho=float(cumulative),
              heavy_L1_above_cut_energy_over_rho=float(Emeasure@abs(lag)),
              heavy_stationary_above_cut_energy_over_rho=float(Emeasure@b[0])))
        print(order,pieces,j["time"],file=sys.stderr,flush=True)
    return dict(scope="Driven high-parent correction on a specified piecewise entropy path, not feedback or a bound.",
                order=order,interpolation_pieces_per_archived_interval=pieces,
                parent_initial_condition=initial,cut_GeV=C,rho_GeV4=float(rho),H_GeV=float(H),
                source_hashes=[x["hash"] for x in raw],
                initial_prompt=float(rate_factor@arr[0][3]),rows=moments)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("output",type=Path)
    p.add_argument("states",type=Path,nargs="+")
    p.add_argument("--order",type=int,default=6)
    p.add_argument("--pieces",type=int,default=1)
    p.add_argument("--initial",choices=("equilibrium","vacuum"),default="equilibrium")
    a=p.parse_args()
    assert a.order>=2 and 1<=a.pieces<=8 and len(a.states)>1
    atomic(a.output,run(a.eos,a.states,a.order,a.pieces,a.initial))
