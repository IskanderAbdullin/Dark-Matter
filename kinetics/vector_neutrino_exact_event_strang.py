"""Exact single-reaction Riccati flow and symmetric Strang quantum event splitting.

For an on-shell event Zprime_j <-> nu_i+antinu_l:
J=b(1-f)(1-a)-(1+b)fa = b(1-f-a)-fa.
With event extent x, b=B-x/wV, f=F+x/wF, a=A+x/wA:
J(x)=q0+q1*x+q2*x² exactly. Each event therefore solves
dx/dt=W*J(x) analytically, with bounded Fermi/Bose occupations,
physical C/L/E and nondecreasing entropy. Palindromic forward
half-events followed by backward-order half-events is symmetric
second-order Strang for this FINITE graph (not physical Γelastic).
Stiff regimes still require small dt/cutoff and independent tests.
"""
from __future__ import annotations
import argparse
import json
import math
import numpy as np
from kinetics.vector_neutrino_shell_radial_endpoints import adjusted_graph
from kinetics.vector_neutrino_reversible_event_graph import (
    synthetic_spectra,ledger,
)


def exact_extent(B,F,A,wV,wF,wA,W,dt):
    """Integrate exact local quantum gain-loss Riccati from x(0)=0.

    All weights must be positive; the explicit root is evaluated
    with expm1 for dt near zero. This is not implicit-Euler.
    """
    if not (all(math.isfinite(z) for z in (B,F,A,wV,wF,wA,W,dt))
            and B>=0 and 0<=F<=1 and 0<=A<=1
            and min(wV,wF,wA,W)>0 and dt>=0):
        raise ValueError("positive weights and physical Bose/Fermi occupations required")
    if dt==0:return 0.
    q0=B*(1-F-A)-F*A
    if q0==0:return 0.
    q1=-(B*(1/wF+1/wA)+(1-F-A)/wV+F/wA+A/wF)
    q2=(1/wF+1/wA)/wV-1/(wF*wA)
    discr=q1*q1-4*q2*q0
    if discr<0:
        if discr> -1e-13*max(q1*q1,1e-300):
            discr=0.
        else:
            raise ArithmeticError("negative Riccati discriminant")
    D=math.sqrt(discr)
    denom=D-q1
    if not(D>0 and abs(denom)>1e-300):
        raise ArithmeticError("degenerate exact event roots: requires independent solver")
    stable=2*q0/denom
    ratio=(-q1-D)/denom
    decay=math.exp(-D*W*dt) if D*W*dt<740 else 0.
    delta=-math.expm1(-D*W*dt) if D*W*dt<740 else 1.
    x=stable*delta/(1-ratio*decay)
    lo=max(-F*wF,-A*wA)
    hi=min(B*wV,(1-F)*wF,(1-A)*wA)
    tol=2e-13*max(abs(hi),abs(lo),1e-15)
    if not (lo-tol<=x<=hi+tol):
        raise ArithmeticError("analytic event would leave physical occupation domain")
    return x


def single_event(g,b,fn,fa,event,dt):
    v,i,l,W=event
    wV,wF=g["wV"],g["wnu"]
    x=exact_extent(b[v],fn[i],fa[l],wV[v],wF[i],wF[l],W,dt)
    b[v]-=x/wV[v]
    fn[i]+=x/wF[i]
    fa[l]+=x/wF[l]


def exact_sweep(g,b,fn,fa,dt,*,scheme="strang"):
    if scheme not in ("strang","forward","reverse"):
        raise ValueError("scheme")
    b=np.array(b,float,copy=True)
    fn=np.array(fn,float,copy=True)
    fa=np.array(fa,float,copy=True)
    if (b.shape!=g["wV"].shape or fn.shape!=g["wnu"].shape or
        fa.shape!=g["wnu"].shape or not np.all(np.isfinite(b)) or
        not np.all(np.isfinite(fn)) or not np.all(np.isfinite(fa)) or
        np.min(b)<0 or np.min(fn)<0 or np.min(fa)<0 or
        np.max(fn)>1 or np.max(fa)>1):
        raise ValueError("physical nonnegative Bose and bounded Pauli distributions required")
    edges=g["edges"]
    if scheme=="strang":
        for ev in edges:single_event(g,b,fn,fa,ev,dt/2)
        for ev in reversed(edges):single_event(g,b,fn,fa,ev,dt/2)
    elif scheme=="forward":
        for ev in edges:single_event(g,b,fn,fa,ev,dt)
    else:
        for ev in reversed(edges):single_event(g,b,fn,fa,ev,dt)
    return b,fn,fa


def pilot(*,h=.2,m=.4,emax=4.,total_tau=4.,dt=.0625,
          scheme="strang"):
    if not (total_tau>=0 and dt>0 and all(math.isfinite(x) for x in
             (total_tau,dt)) and
            abs(round(total_tau/dt)*dt-total_tau)<1e-10):
        raise ValueError("time and dt must be finite and commensurate")
    g=adjusted_graph(m,h,emax,scheme="shell_cell_integral")
    initial=synthetic_spectra(g)
    before=ledger(g,*initial)
    b,fn,fa=(x.copy() for x in initial)
    min_step_dS=math.inf
    last=before
    for _ in range(round(total_tau/dt)):
        b,fn,fa=exact_sweep(g,b,fn,fa,dt,scheme=scheme)
        state=ledger(g,b,fn,fa)
        min_step_dS=min(min_step_dS,state["entropy"]-last["entropy"])
        last=state
    return (b,fn,fa),{
      "h":h,"total_tau":total_tau,"dt":dt,
      "scheme":scheme,"steps":round(total_tau/dt),
      "low_vector_b":float(b[0]),
      "total_entropy_gain":last["entropy"]-before["entropy"],
      "min_step_entropy_gain":min_step_dS,
      "max_C_drift":abs(last["reduced_charge"]-before["reduced_charge"]),
      "max_L_drift":abs(last["lepton_asymmetry"]-before["lepton_asymmetry"]),
      "max_E_drift":abs(last["energy"]-before["energy"]),
      "min_b":float(np.min(b)),
      "min_fermion":float(min(np.min(fn),np.min(fa))),
      "max_fermion":float(max(np.max(fn),np.max(fa))),
      "synthetic_only_no_cosmological_gamma_or_86eV_DM":True,
    }


if __name__=="__main__":
    a=argparse.ArgumentParser()
    a.add_argument("--h",type=float,default=.2)
    a.add_argument("--tau",type=float,default=4.)
    a.add_argument("--dt",type=float,default=.0625)
    p=a.parse_args()
    rows=[]
    reference,_=pilot(h=p.h,total_tau=p.tau,dt=p.dt/8)
    for dt in (p.dt*8,p.dt*4,p.dt*2,p.dt):
        if abs(round(p.tau/dt)*dt-p.tau)>1e-10:continue
        r,info=pilot(h=p.h,total_tau=p.tau,dt=dt)
        info["maximum_spectral_abs_error_vs_ref"]=max(
            float(np.max(np.abs(a-b))) for a,b in zip(r,reference))
        rows.append(info)
    print(json.dumps({"reference_dt":p.dt/8,"scan":rows,
        "not_a_physical_continuum_relaxation_rate":True},indent=2))
