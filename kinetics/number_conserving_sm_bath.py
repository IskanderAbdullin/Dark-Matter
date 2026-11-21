"""Energy-accounted, number-conserving thermal elastic bath for massive FD chi1.

This is a phenomenological relaxation-time (BGK) operator, NOT a derived
microphysical Standard Model scattering amplitude or rethermalization rate.

On a fixed COMOVING momentum grid q=a p and at a fixed instantaneous
scale factor a, the bath tends to physical FD:
   f_eq(q) = 1/(exp[(sqrt((q/a)^2+M²)-mu)/T_bath]+1).
Crucially, chemical mu is determined at EVERY time from the actual
mother number Nchi1 = g/(2pi²) sum[w(q) q² f(q)] so that the target
has exactly the SAME comoving number. This models *elastic* scattering,
not creation/annihilation of chi1 by an untracked reservoir.

Collision operator: C_el[f] = gamma_el*(f_eq[f]-f).
It is number-conserving in the same finite quadrature on which it
acts. At frozen a,Tbath, the exact substep is the convex mixture
  f_after = exp(-gamma*dt)*f_before
            + (1-exp(-gamma*dt))*f_eq.
This guarantees 0<=f<=1 and exact preservation of Nchi1 on each
substep; the relative Fermi free-energy/KL functional decreases.

Every nonzero dark-sector *physical* energy transfer is tracked as
   Delta E_com,dark = g/(2pi²) sum[w q² sqrt(q²+(a M)²) Delta f]
   Delta E_com,SM = -Delta E_com,dark.
The SM bath temperature is NOT yet evolved from this transfer.
For evolving a(t), pressure work and radiation entropy have to be
tracked separately, not hidden in this zero-sum energy bookkeeping.

The pinned BEST solver supports adding (gain,loss) custom collisions,
but using it does not by itself implement a physical SM heat reservoir.
This module provides an independently testable conservative building
block before coupling Q to a dynamic radiation bath.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import expit, xlogy


def _weights_and_state(q,f,g):
    q=np.asarray(q,dtype=float)
    f=np.asarray(f,dtype=float)
    if (q.ndim!=1 or f.ndim!=1 or len(q)<4 or len(q)!=len(f)
            or not np.all(np.isfinite(q)) or not np.all(np.isfinite(f))
            or np.any(q<=0) or np.any(np.diff(q)<=0)
            or np.any(f<0) or np.any(f>1)
            or not math.isfinite(g) or g<=0):
        raise ValueError("Require increasing positive finite q, fermion occupations 0..1, and g>0")
    weights=np.zeros_like(q)
    widths=np.diff(q)/2
    weights[:-1]+=widths
    weights[1:]+=widths
    return q,f,g/(2*np.pi**2)*weights*q*q


def _projected_target(q,f,a,mass,Tbath,g):
    q,f,W=_weights_and_state(q,f,g)
    if (not all(math.isfinite(x) for x in (a,mass,Tbath))
            or a<=0 or mass<=0 or Tbath<=0):
        raise ValueError("Require positive finite a, fermion mass, and thermal bath T")
    physical_E=np.hypot(q/a,mass)
    current_number=float(np.dot(W,f))
    max_number=float(W.sum())
    if current_number<=0:
        return q,f,W,physical_E,np.zeros_like(f),-math.inf
    if current_number>=max_number:
        return q,f,W,physical_E,np.ones_like(f),math.inf
    def mismatch(mu):
        return float(np.dot(W,expit((mu-physical_E)/Tbath))-current_number)
    mu=brentq(mismatch,float(physical_E.min()-100*Tbath),
                       float(physical_E.max()+100*Tbath),
                       xtol=1e-13,rtol=1e-12)
    eq=expit((mu-physical_E)/Tbath)
    return q,f,W,physical_E,eq,float(mu)


def number_conserving_bath_rates(q,f,a=1.18,mass=.1,Tbath=.82,g=4,gamma=1.2):
    q,f,W,Ephys,eq,mu=_projected_target(q,f,a,mass,Tbath,g)
    if not math.isfinite(gamma) or gamma<0:
        raise ValueError("Require nonnegative finite relaxation rate")
    gain=gamma*eq
    loss=gamma*f
    return {
        "gain":gain,"loss":loss,
        "rate":gain-loss,
        "chemical_mu_number_projected":mu,
        "comoving_number_derivative":float(np.dot(W,gain-loss)),
        "comoving_dark_energy_derivative":float(
            np.dot(W*a*Ephys,gain-loss)),
        "comoving_SM_bath_energy_derivative":float(
            -np.dot(W*a*Ephys,gain-loss)),
        "phenomenological_elastic_BGK_only":True
    }


def _relative_FD_entropy(f,eq,weights):
    # Clipping only protects the LOG of a probability; the actual
    # charge/energy quadratures and particle update are NOT clipped.
    x=np.clip(f,1e-15,1-1e-15)
    y=np.clip(eq,1e-15,1-1e-15)
    return float(np.dot(weights,
        xlogy(x,x/y)+xlogy(1-x,(1-x)/(1-y))))


def exact_frozen_bath_step(q,f,a=1.18,mass=.1,Tbath=.82,
                           g=4,gamma=1.2,dt=.4):
    q,f,W,Ephys,eq,mu=_projected_target(q,f,a,mass,Tbath,g)
    if not all(math.isfinite(x) for x in (gamma,dt)) or gamma<0 or dt<0:
        raise ValueError("Require nonnegative finite collision rate and step")
    fraction=-math.expm1(-gamma*dt)
    after=f+fraction*(eq-f)
    before_num=float(np.dot(W,f))
    after_num=float(np.dot(W,after))
    before_energy=float(np.dot(W*a*Ephys,f))
    after_energy=float(np.dot(W*a*Ephys,after))
    delta_dark=after_energy-before_energy
    return {
        "q":q,"before":f,"after":after,"target_fermi_dirac":eq,
        "bath_T":Tbath,"chemical_mu_number_projected":mu,
        "scale_factor":a,"mass":mass,"gamma_relax":gamma,"dt":dt,
        "number_before":before_num,"number_after":after_num,
        "comoving_number_residual":after_num-before_num,
        "dark_comoving_energy_change":delta_dark,
        "SM_bath_comoving_energy_change":-delta_dark,
        "dark_plus_bath_comoving_energy_exchange_residual":0.,
        "fermionic_KL_before":_relative_FD_entropy(f,eq,W),
        "fermionic_KL_after":_relative_FD_entropy(after,eq,W),
        "minimum_occupation":float(after.min()),
        "maximum_occupation":float(after.max()),
        "elastic_chi1_number_conserved":True,
        "bath_T_not_self_consistently_evolved":True,
        "not_a_derived_SM_2to2_matrix_element":True,
        "not_a_converged_cosmological_DM_model":True
    }


def demonstration():
    q=np.geomspace(1e-4,30.,96)
    f=np.clip(expit(-np.hypot(q,.1))*
              (1+.05*np.tanh(np.log1p(q))),0,1)
    outcome=exact_frozen_bath_step(q,f)
    return {key:value for key,value in outcome.items()
            if not isinstance(value,np.ndarray)}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",default="outputs/number_conserving_SM_BGK_ledger.json")
    args=parser.parse_args()
    r=demonstration()
    path=Path(args.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(r,indent=2),encoding="utf-8")
    print(json.dumps(r,indent=2))
