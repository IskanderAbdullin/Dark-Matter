"""Feasibility and optional bounded two-moment correction for massive parent grids.

At a given x=M/T, let N_i=C_i f_i, where C_i is the q²dq particle
capacity of an energy shell and 0<=f_i<=1. The shell mean physical
energy is E_i=<sqrt(q²+x²)>_{q²dq}. E_i strictly increases with q.

At fixed total N, allowed energy is exactly between filling LOWEST
energy cells first and filling HIGHEST energy cells first (rearrangement
inequality). This is an independent feasibility bound, not a model
of quantum collisions.

Optionally a number-conservative Pauli-preserving remap can repair its
finite-volume energy defect by transferring minimal number between low
occupied and high unfilled cells (or reverse). This global artificial
redistribution is NOT collisionless Liouville dynamics: it changes f(q),
may change entropy, and must never silently be included in a physical
FRW collision solver or used as CLASS initial data.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre

from parent_q_remap import edges_checked, rebin, number, energy, entropy


def shell_capacities_and_energy(q_edges, x, g=4, gauss_n=24):
    q=edges_checked(q_edges)
    if x<=0 or g<=0 or gauss_n<3:
        raise ValueError("Positive x, g, Gauss order required")
    capacity=g/(6*np.pi**2)*np.diff(q**3)
    z,w=roots_legendre(gauss_n)
    low=q[:-1,None];high=q[1:,None]
    sample=.5*(low+high)+.5*(high-low)*z
    weighted=.5*(high-low)*w*sample**2
    eavg=(weighted*np.hypot(sample,x)).sum(axis=1)/(np.diff(q**3)/3)
    if np.any(np.diff(eavg)<=0):
        raise ArithmeticError("Shell mean energy must increase")
    return capacity,eavg


def _extremal_energy(capacity,ebar,particles,reverse=False):
    remaining=particles
    E=0.
    order=range(len(capacity)-1,-1,-1) if reverse else range(len(capacity))
    for i in order:
        allocated=min(capacity[i],remaining)
        E+=allocated*ebar[i]
        remaining-=allocated
        if remaining<=1e-13*max(1,particles):
            break
    if remaining>1e-11*max(1,particles):
        raise ArithmeticError("Could not fill requested particles")
    return float(E)


def feasible_energy_interval(q_edges,x,particles,g=4):
    capacity,ebar=shell_capacities_and_energy(q_edges,x,g)
    if particles<0 or particles>capacity.sum()*(1+1e-12):
        raise ValueError("Requested particle number exceeds Pauli capacity")
    particles=min(float(particles),float(capacity.sum()))
    return {
        "minimum_energy":_extremal_energy(capacity,ebar,particles),
        "maximum_energy":_extremal_energy(capacity,ebar,particles,reverse=True),
        "particle_capacity":float(capacity.sum()),
        "requested_particles":particles,
        "statement":"Bounds follow from occupancy 0<=f<=1 and strictly increasing shell energies"
    }


def repair_energy_by_bounded_transfer(old_edges,new_edges,old_f,x,tolerance=2e-11):
    """Optional artificial correction; exact N/E, not exact Liouville shape."""
    new=edges_checked(new_edges)
    f0=rebin(old_edges,new,old_f)
    capacity,ebar=shell_capacities_and_energy(new,x)
    N=capacity*f0
    old_energy=energy(old_edges,old_f,x,gauss_n=24)
    initial_new_energy=float(np.dot(N,ebar))
    old_particles=number(old_edges,old_f)
    interval=feasible_energy_interval(new,x,old_particles)
    scale=max(1.,old_energy)
    if not (interval["minimum_energy"]-tolerance*scale <= old_energy <=
            interval["maximum_energy"]+tolerance*scale):
        raise ValueError("Exact energy at fixed number infeasible under Pauli bounds")
    deficit=old_energy-initial_new_energy
    i,j=0,len(N)-1
    transferred=0.
    while abs(deficit)>tolerance*scale:
        if deficit>0:
            while i<j and N[i]<=1e-16:i+=1
            while j>i and capacity[j]-N[j]<=1e-16:j-=1
            if i>=j:break
            unit=ebar[j]-ebar[i]
            step=min(N[i],capacity[j]-N[j],deficit/unit)
            if step<=0:break
            N[i]-=step;N[j]+=step
            deficit-=step*unit
        else:
            while j>i and N[j]<=1e-16:j-=1
            while i<j and capacity[i]-N[i]<=1e-16:i+=1
            if i>=j:break
            unit=ebar[j]-ebar[i]
            step=min(N[j],capacity[i]-N[i],-deficit/unit)
            if step<=0:break
            N[j]-=step;N[i]+=step
            deficit+=step*unit
        transferred+=step
    corrected=np.clip(N/capacity,0.,1.)
    E2=float(np.dot(capacity*corrected,ebar))
    particle_error=number(new,corrected)-old_particles
    if abs(E2-old_energy)>tolerance*scale:
        raise ArithmeticError("Pauli-safe correction cannot match target energy")
    if abs(particle_error)>tolerance*max(1.,old_particles):
        raise ArithmeticError("Pauli-safe correction lost particle number")
    before_entropy=entropy(old_edges,old_f)
    after_entropy=entropy(new,corrected)
    return {
        "corrected_occupations":corrected,
        "number_error":particle_error,
        "relative_energy_error":E2/old_energy-1,
        "raw_remap_relative_energy_error":initial_new_energy/old_energy-1,
        "shifted_particles_fraction":transferred/max(old_particles,1e-300),
        "numeric_entropy_change":after_entropy-before_entropy,
        "feasible_energy_interval":interval,
        "classification":"BOUNDED FEASIBILITY CONSTRUCTION ONLY, NOT PHYSICAL LIOUVILLE"
    }


def run():
    old=np.r_[0.,np.geomspace(.01,10.,100)]
    new=np.r_[0.,np.geomspace(.012,10.,101)]
    f=1/(1+np.exp(np.hypot(.5*(old[1:]+old[:-1]),.1)))
    x=.2
    report=repair_energy_by_bounded_transfer(old,new,f,x)
    return {k:(list(map(float,v)) if k=="corrected_occupations" else v)
            for k,v in report.items()}


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="outputs/pauli_energy_feasibility.json")
    o=p.parse_args()
    result=run()
    dst=Path(o.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,indent=2),encoding="utf-8")
    print(json.dumps(result,indent=2))
