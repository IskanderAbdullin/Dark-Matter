"""Independent constrained-entropy FD/BE equilibrium for the corrected finite lattice.

Z' <-> nu + antinu, EXACT three physical collisional invariants on the
half-offset daughter grid at tested (m,h,Emax); all time evolution
uses the existing positivity-preserving microreversible event solver
without resetting either neutrino species to FD.

This is a finite-dimensional NUMERICAL TEST, not early-universe
thermalization, a continuous strong equation or physical chi gamma_el.
"""
from __future__ import annotations
import argparse
import json
import math
import numpy as np
from scipy.optimize import root
from scipy.special import expit
from kinetics.vector_neutrino_shell_radial_endpoints import adjusted_graph
from kinetics.vector_neutrino_reversible_event_graph import (
    synthetic_spectra,ledger,sweep,
)


def constrained_equilibrium(g,initial,*,initial_guess=(0.,-.1,-.1)):
    """Maximum-entropy occupation constrained by C,L,E, with ν asymmetry.

    Unknowns x=(log T, μν/T, μantiν/T) remove positivity of T.
    Physical Bose state demands E_j/T−x1−x2>0 for all vector modes.
    This independent root does not call the collision or time integrator.
    """
    C=ledger(g,*initial)
    target=np.array([C["reduced_charge"],C["lepton_asymmetry"],C["energy"]])
    E,k,wv,wn=(g[a] for a in ("E","k","wV","wnu"))

    def model(x,with_state=False):
        if not np.all(np.isfinite(x)) or abs(x[0])>15:
            return None if with_state else np.full(3,1e10)
        invt=math.exp(-x[0])
        argument=E*invt-x[1]-x[2]
        if min(argument)<=0:
            return None if with_state else np.full(3,1e10)
        b=1/np.expm1(argument)
        f=expit(x[1]-k*invt)
        a=expit(x[2]-k*invt)
        if with_state:return b,f,a
        result=np.array([np.dot(wn,f+a)+2*np.dot(wv,b),
                         np.dot(wn,f-a),
                         np.dot(wn,k*(f+a))+np.dot(wv,E*b)])
        return (result-target)/np.array([max(abs(target[0]),.1),.1,
                                          max(abs(target[2]),.1)])
    res=root(model,np.asarray(initial_guess,float),method="hybr",tol=2e-12)
    if not res.success and np.max(np.abs(res.fun))>5e-12:
        raise RuntimeError("constrained-equilibrium root not numerically certified")
    final=model(res.x,True)
    if final is None:raise RuntimeError("invalid Bose equilibrium state")
    check=ledger(g,*final)
    residual=max(abs(check["reduced_charge"]-target[0]),
                 abs(check["lepton_asymmetry"]-target[1]),
                 abs(check["energy"]-target[2]))
    if residual>2e-11:raise ArithmeticError("conserved quantities do not match")
    return {"T":float(math.exp(res.x[0])),
            "mu_nu_over_T":float(res.x[1]),
            "mu_antinu_over_T":float(res.x[2]),
            "muV_over_T":float(res.x[1]+res.x[2]),
            "residual_max":float(residual),
            "entropy_start":float(C["entropy"]),
            "entropy_equilibrium":float(check["entropy"]),
            "entropy_gain_max":float(check["entropy"]-C["entropy"]),
            "nonzero_lepton_asymmetry":abs(target[1])>0,
            "equilibrium_is_discrete_and_not_actual_SM":True},final


def long_time_pilot(*,h=.2,m=.4,emax=4.,tau=500.,dt=1.):
    if not (tau>=0 and dt>0 and abs(round(tau/dt)*dt-tau)<1e-9):
        raise ValueError("finite time interval must divide into dt")
    g=adjusted_graph(m,h,emax,scheme="shell_cell_integral")
    initial=synthetic_spectra(g)
    eq,ref=constrained_equilibrium(g,initial)
    b,f,a=(v.copy() for v in initial)
    start=ledger(g,b,f,a)
    minS=math.inf
    previousS=start["entropy"]
    for _ in range(round(tau/dt)):
        b,f,a=sweep(g,b,f,a,dt)
        now=ledger(g,b,f,a)["entropy"]
        minS=min(minS,now-previousS)
        previousS=now
    end=ledger(g,b,f,a)
    error=max(float(np.max(np.abs(x-y))) for x,y in
              zip((b,f,a),ref))
    return {"h":h,"m":m,"emax":emax,"tau_Gamma0":tau,
            "dt_Gamma0":dt,"steps":round(tau/dt),
            "microreversible_edges":len(g["edges"]),
            "conserved_C_residual":end["reduced_charge"]-start["reduced_charge"],
            "conserved_L_residual":end["lepton_asymmetry"]-start["lepton_asymmetry"],
            "conserved_E_residual":end["energy"]-start["energy"],
            "entropy_gain":end["entropy"]-start["entropy"],
            "entropy_gap_to_independent_constrained_max":
                eq["entropy_equilibrium"]-end["entropy"],
            "min_step_entropy_gain":minS,
            "max_occupation_abs_error_to_FD_BE_equilibrium":error,
            "lowest_vector_b":float(b[0]),
            "nu_and_antinu_not_reset_to_thermal":True,
            "physical_gamma_el_or_86eV_relic_validated":False}


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--h",type=float,default=.2)
    p.add_argument("--tau",type=float,default=500.)
    p.add_argument("--dt",type=float,default=1.)
    a=p.parse_args()
    g=adjusted_graph(.4,a.h,4.,scheme="shell_cell_integral")
    eq,_=constrained_equilibrium(g,synthetic_spectra(g))
    print(json.dumps({"equilibrium":eq,"time":long_time_pilot(
        h=a.h,tau=a.tau,dt=a.dt)},indent=2))
