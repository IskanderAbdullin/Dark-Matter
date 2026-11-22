"""Exact angular CELL OVERLAP for energy-matched Zprime->nu+antinu graph.

Complement to reversible_event_graph: midpoint neutrino energies
k_i=(i+1/2)h and vector E_j=jh, so E_j=k_i+k_(j-i-1) exactly.
For each vector shell the physical k range [(E-p)/2,(E+p)/2]
overlaps ν cell [ih,(i+1)h] by a positive analytic length
delta_k; sum(delta_k)=p exactly. In contrast, center-hit weights
over/undercount even a constant angular integrand.
Still midpoint energy/finite-Emax approximation, not continuum solver.
"""
from __future__ import annotations
import math
import numpy as np
from kinetics.vector_neutrino_reversible_event_graph import (
    synthetic_spectra,event_moments,independent_continuum_weak,
)

def graph_overlap(m=.4,h=.2,emax=4.,gamma0=1.,*,cell_overlap=True):
    if not(all(math.isfinite(v) and v>0 for v in (m,h,emax,gamma0))
           and emax>m and abs(emax/h-round(emax/h))<1e-10):
        raise ValueError("positive scales and Emax/h integer")
    jvec=np.arange(int(math.floor(m/h))+1,int(round(emax/h))+1,dtype=int)
    E=jvec*h
    p=np.sqrt(E*E-m*m)
    k=(np.arange(int(round(emax/h)),dtype=float)+.5)*h
    wV=3/(2*math.pi**2)*E*p*h
    wnu=2/(2*math.pi**2)*k*k*h
    edges=[]
    for v,(j,pp) in enumerate(zip(jvec,p)):
        klo=(E[v]-pp)/2
        khi=(E[v]+pp)/2
        for i in range(int(j)):
            if cell_overlap:
                d=max(0.,min((i+1)*h,khi)-max(i*h,klo))
            else:
                center=(i+.5)*h
                d=h if klo<=center<=khi else 0.
            if d>0:
                W=wV[v]*gamma0*m/E[v]*d/pp
                edges.append((v,i,int(j-i-1),W))
    return {"E":E,"p":p,"k":k,"wV":wV,"wnu":wnu,
            "edges":edges,"m":m,"h":h,"emax":emax,
            "gamma0":gamma0,"cell_overlap":cell_overlap}

def shell_angular_normalization(g):
    """For flat quantum event factor, total decay polar probability =1."""
    sums=np.zeros(len(g["E"]))
    for v,i,l,W in g["edges"]:
        sums[v]+=W
    exact=g["wV"]*g["gamma0"]*g["m"]/g["E"]
    return sums/exact

if __name__=="__main__":
    import json
    ref=independent_continuum_weak(.4,4)
    rows=[]
    for h in (.4,.2,.1,.05):
        record={"h":h}
        for label,overlap in (("center",False),("cell_overlap",True)):
            g=graph_overlap(.4,h,4,cell_overlap=overlap)
            x=event_moments(g,*synthetic_spectra(g))
            record[label]={
                "edges":len(g["edges"]),
                "max_constant_angular_normalization_error":
                    float(np.max(abs(shell_angular_normalization(g)-1))),
                "dNV_relative_to_continuum":x["dNV"]/ref["dNV"]-1,
                "dEV_relative_to_continuum":x["dEV"]/ref["dEV"]-1,
                "entropy_relative_to_continuum":
                    x["entropy_production"]/ref["entropy_production"]-1}
        rows.append(record)
    print(json.dumps({"reference":ref,"grid_scan":rows,
                      "no_continuum_accuracy_certificate":True},indent=2))
