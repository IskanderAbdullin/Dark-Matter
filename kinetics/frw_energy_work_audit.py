"""Physical FRW pressure-work integral and separate numerical energy ledger.

For freely streaming massive f(q) at constant gstarS and fixed q=p/T,
the dimensionless physical mother energy rho_bar(x)=g/(2pi²)*integral
q² f(q)*sqrt(q²+x²)dq obeys

 d rho_bar/dx = g/(2pi²)*integral q² f(q)*x/E dq
              = (rho_bar-3P_bar)/x.

The integral of this differential work over x exactly equals the
change of rho_bar for an unchanged phase-space distribution.

This independent Gauss-in-(q,x) check identifies legitimate FRW pressure
work separately from (a) finite-cell collision physical-energy
mismatch and (b) artificial number-conserving remap energy drift.
A bookkeeping ledger that closes numerically does NOT prove that
(a)+(b) vanish, or that the dark sector is physically energy-closed.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre, expit

from parent_q_remap import edges_checked, energy, thermal_shell_averages


def physical_parent_pressure(q_edges,f,x,g=4,gauss_n=24):
    q=edges_checked(q_edges)
    f=np.asarray(f,dtype=float)
    if x<=0 or f.shape!=(len(q)-1,) or g<=0:
        raise ValueError("Invalid physical pressure inputs")
    nodes,w=roots_legendre(gauss_n)
    low=q[:-1,None];high=q[1:,None]
    qnodes=.5*(low+high)+.5*(high-low)*nodes
    base=.5*(high-low)*w*qnodes**2*f[:,None]
    E=np.hypot(qnodes,x)
    rho=float(g/(2*np.pi**2)*np.sum(base*E))
    pressure=float(g/(6*np.pi**2)*np.sum(base*qnodes**2/E))
    derivative=float(g/(2*np.pi**2)*np.sum(base*x/E))
    return {
        "rho_over_T4":rho,
        "pressure_over_T4":pressure,
        "d_rho_over_T4_dx":derivative,
        "pressure_continuity_residual":x*derivative-rho+3*pressure
    }


def direct_frw_pressure_work(q_edges,f,x0,x1,g=4,q_order=24,x_order=64):
    """Independent x-integral of fixed-q parent mass work, no collisions."""
    q=edges_checked(q_edges)
    f=np.asarray(f,dtype=float)
    if f.shape!=(len(q)-1,) or x0<=0 or x1<x0:
        raise ValueError("Invalid expansion interval/distribution")
    if x0==x1:return 0.
    nq,wq=roots_legendre(q_order)
    nx,wx=roots_legendre(x_order)
    left=q[:-1,None];right=q[1:,None]
    Q=.5*(left+right)+.5*(right-left)*nq
    measure=.5*(right-left)*wq*Q**2*f[:,None]
    X=.5*(x0+x1)+.5*(x1-x0)*nx
    # Integrate x/E(q,x) with Gauss-Legendre in x INDEPENDENTLY of
    # evaluating the two endpoint physical energies.
    integrated=np.sum(wx[None,None,:]*X[None,None,:]/
                      np.hypot(Q[:,:,None],X[None,None,:]),axis=-1)
    return float(g/(2*np.pi**2)*.5*(x1-x0)*np.sum(measure*integrated))


def verify_free_expansion(q_edges,f,x0,x1):
    work=direct_frw_pressure_work(q_edges,f,x0,x1)
    difference=energy(q_edges,f,x1,gauss_n=24)-energy(q_edges,f,x0,gauss_n=24)
    residual=work-difference
    if abs(residual)>3e-10*max(1,abs(difference)):
        raise ArithmeticError("Independent FRW work integral disagrees with physical energy change")
    return {
        "x0":x0,"x1":x1,
        "integral_of_rho_minus_3P_over_x":work,
        "direct_physical_parent_energy_change":difference,
        "independent_integration_residual":residual,
        "physical_status":"FREE MASSIVE PARENT IN PRESCRIBED CONSTANT-GSTAR FRW"
    }


def audit_split_ledger(simulation):
    steps=simulation["steps"]
    if not steps:
        raise ValueError("Missing split-step evidence")
    physical=sum(float(x["physical_energy_change_from_start_to_end"]) for x in steps)
    work=sum(float(x.get("physical_free_expansion_energy_change",
                         x.get("collisionless_FRW_energy_change"))) for x in steps)
    errors=sum(float(x["collision_plus_remap_non_FRW_energy_error"]) for x in steps)
    closure=physical-work-errors
    maximum_individual=max(abs(float(x["physical_energy_ledger_residual"])) for x in steps)
    if abs(closure)>1e-10 or maximum_individual>1e-10:
        raise ArithmeticError("Energy ledger decomposition did not close")
    return {
        "physical_energy_change_in_T4_units":physical,
        "legitimate_collisionless_expansion_work_in_T4_units":work,
        "collision_shell_and_remap_energy_defects_in_T4_units":errors,
        "ledger_residual":closure,
        "maximum_individual_ledger_residual":maximum_individual,
        "unphysical_energy_error_is_not_assumed_zero":True,
        "not_a_self_consistent_cosmology":True
    }


def run():
    q=np.r_[0.,np.geomspace(.01,30.,120)]
    f=thermal_shell_averages(q,.05)
    rows=[verify_free_expansion(q,f,.05,x) for x in (.1,.5,2.,10.)]
    return {
        "scope":"Analytic FRW Liouville pressure work vs independent x quadrature",
        "model":"collisionless massive fermion, fixed comoving q and constant gstarS",
        "cases":rows,
        "warning":"The independently tested pressure identity is an essential component but not proof of full collision+FRW energy conservation, SM energy closure or Bose dark matter."
    }


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="outputs/frw_energy_work_audit.json")
    opts=p.parse_args()
    result=run()
    dest=Path(opts.output)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,indent=2),encoding="utf-8")
    print(json.dumps(result,indent=2))
