"""Vector-energy endpoint quadrature for energy-matched Zprime->nu+antinu events.

The exact daughter angular cell overlap removes angular alias, but its
original wV(E_j)=gV p_j E_j h/(2pi²) assigns FULL h to E_j=Emax,
while independent continuum weak reference integrates E<=Emax.
Two physically motivated finite-interval variants:
- upper_trapezoid: half-weight at the upper endpoint;
- shell_cell_integral: exact integral of E sqrt(E²-mV²) within
  intersected E-bin, with E_j retained as discrete collocation energy.
Both preserve positive events and exact DISCRETE energy and H theorem.
Neither corrects missing threshold bin or actual continuum energy
moment collocation, nor is a physical cosmological solver.
"""
from __future__ import annotations
import math
import numpy as np
from kinetics.vector_neutrino_shell_cell_overlap import graph_overlap,shell_angular_normalization
from kinetics.vector_neutrino_reversible_event_graph import (
    event_moments,synthetic_spectra,independent_continuum_weak,
)


def adjusted_graph(m=.4,h=.2,emax=4.,*,scheme="upper_trapezoid"):
    if scheme not in ("original","upper_trapezoid","shell_cell_integral"):
        raise ValueError("unknown physical radial integration scheme")
    g=graph_overlap(m,h,emax,cell_overlap=True)
    E=g["E"];p=g["p"];old=g["wV"]
    fac=np.ones_like(old)
    if scheme=="upper_trapezoid":
        fac[-1]=.5
    elif scheme=="shell_cell_integral":
        lo=np.maximum(m,E-h/2)
        hi=np.minimum(emax,E+h/2)
        wnew=3/(2*math.pi**2)/3*(
            np.maximum(0,hi*hi-m*m)**1.5-
            np.maximum(0,lo*lo-m*m)**1.5)
        fac=wnew/old
    g["wV"]=old*fac
    g["edges"]=[(v,i,l,W*fac[v]) for v,i,l,W in g["edges"]]
    g["vector_radial_scheme"]=scheme
    return g


def scan(m=.4,emax=4.,steps=(.4,.2,.1,.05)):
    ref=independent_continuum_weak(m,emax)
    out=[]
    for h in steps:
        row={"h":h}
        for label in ("original","upper_trapezoid","shell_cell_integral"):
            g=adjusted_graph(m,h,emax,scheme=label)
            x=event_moments(g,*synthetic_spectra(g))
            row[label]={
                "dNV_rel":x["dNV"]/ref["dNV"]-1,
                "dEV_rel":x["dEV"]/ref["dEV"]-1,
                "entropy_rel":x["entropy_production"]/ref["entropy_production"]-1,
                "flat_angular_shell_max_abs_err":
                  float(np.max(abs(shell_angular_normalization(g)-1)))}
        out.append(row)
    return {"continuum":ref,"scan":out,
            "physical_continuum_time_not_verified":True}


if __name__=="__main__":
    import json
    print(json.dumps(scan(),indent=2))
