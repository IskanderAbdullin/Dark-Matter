"""Quantum Jensen entropy audit for conservative Bose/Pauli UV remapping.

The quantum per-mode entropies
  s_B(f)=(1+f)log(1+f)-f log f, f>=0
  s_F(f)=-f log f-(1-f)log(1-f), 0<=f<=1
are concave. Exact q³-volume-weighted averaging (with zero-occupation
extension in newly added UV domain) therefore cannot LOWER the
coarse-grained quantum entropy of any one species.

This is a MATHEMATICAL numerical mixing effect. It is NOT physical
dissipation or entropy created by a microscopic collision process.
Conservation of charge and entropy monotonicity are still insufficient
to demonstrate conservation of physical energy and f(q) convergence.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import xlogy

from daughter_uv_regrid import positive_species_regrid,remap_three_species,parent_qedges
from sparse_kinetic_jacobian import SparseClosedDecay
from resonant_graph_repair import bridged_indices
from scipy.special import expit


def quantum_entropy(q_edges,occupation,statistics,g):
    q=np.asarray(q_edges,dtype=float)
    f=np.asarray(occupation,dtype=float)
    if (q.ndim!=1 or len(q)<2 or q[0]!=0 or np.any(np.diff(q)<=0)
            or f.shape!=(len(q)-1,) or not np.all(np.isfinite(f))
            or np.any(f<0) or g<=0):
        raise ValueError("Invalid q-occupation entropy specification")
    if statistics=="bose":
        h=(1+f)*np.log1p(f)-xlogy(f,f)
    elif statistics=="fermi":
        if np.any(f>1):raise ValueError("Pauli violation")
        h=-(xlogy(f,f)+xlogy(1-f,1-f))
    else:
        raise ValueError("Quantum statistics must be bose/fermi")
    return float(g/(6*np.pi**2)*np.dot(np.diff(q**3),h))


def verify_entropy_mixing(old_q,new_q,f,statistics,g):
    S0=quantum_entropy(old_q,f,statistics,g)
    mapped=positive_species_regrid(old_q,new_q,f,statistics,g)
    S1=quantum_entropy(new_q,mapped,statistics,g)
    diff=S1-S0
    if diff < -5e-10*max(1,abs(S0)):
        raise ArithmeticError("Convex remap decreased quantum entropy")
    return {"before":S0,"after":S1,"increase":diff,
            "quantum_statistics":statistics,
            "meaning":"Numerical coarse-graining entropy, NOT a physical collision H-theorem"}


def one_full_repair():
    m0=SparseClosedDecay(bridged_indices(max_index=300),M1_T=.02)
    m1=SparseClosedDecay(bridged_indices(max_index=500),M1_T=.03)
    b=1e4*np.exp(-m0.e/.03)
    f=.2*np.exp(-m0.e)
    parent=expit(-m0.parent_E)
    full=np.r_[b,f,parent]
    oldparent=parent_qedges(m0,.02)
    newparent=parent_qedges(m1,.03)
    species={
       "phi":verify_entropy_mixing(m0.edges,m1.edges,b,"bose",1),
       "chi2":verify_entropy_mixing(m0.edges,m1.edges,f,"fermi",4),
       "chi1":verify_entropy_mixing(oldparent,newparent,parent,"fermi",4),
    }
    remap=remap_three_species(m0,m1,full,.02,.03)
    newstate=remap["remapped_state"]
    b2,f2,p2=m1.unpack(newstate)
    checks={
      "phi":quantum_entropy(m1.edges,b2,"bose",1),
      "chi2":quantum_entropy(m1.edges,f2,"fermi",4),
      "chi1":quantum_entropy(newparent,p2,"fermi",4)
    }
    for key in checks:
        if abs(checks[key]-species[key]["after"])>1e-8:
            raise ArithmeticError("Full 3-species mapping differs from per-species entropy")
    total_diff=sum(d["increase"] for d in species.values())
    if total_diff<-1e-9:
        raise ArithmeticError("Combined quantum entropy decreased")
    return {
        "model":"M1/T=.02 -> .03, daughter index Rmax=300 -> 500, vacuum high-q extension",
        "individual_species_entropy":species,
        "total_quantum_coarse_grained_entropy_increase":float(total_diff),
        "two_charge_conservation_residual":remap["maximum_charge_drift"],
        "pure_regridding_physical_energy_error":
            remap["pure_regridding_physical_energy_error"],
        "numerical_entropy_not_physical_collision_production":True,
        "not_a_validated_cosmological_distribution":True
    }


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="outputs/quantum_regrid_entropy.json")
    opt=p.parse_args()
    data=one_full_repair()
    out=Path(opt.output)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(data,indent=2),encoding="utf-8")
    print(json.dumps(data,indent=2))
