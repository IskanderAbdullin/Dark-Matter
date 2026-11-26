"""All-quartet nonlinear finite χ/e Pauli kinetic evolution at fixed volume.

Rate tensor is derived from full 5D vector-contact scattering but
PAIR-AVERAGED between separately sampled finite-bin event directions.
That finite-lattice closure is NOT the unmodified continuum Boltzmann
operator. Nevertheless all lattice-allowed quartets evolve BOTH χ and
electron occupations without a thermostat, using physically integrated
positive W and a scalar Pauli-bounded implicit reaction step per edge.

Lie-split global evolution preserves separate χ/e numbers, total energy,
Pauli 0<=f<=1 and nondecreasing discrete Fermi entropy, but needs
independent time-order/dt, energy/angular and cutoff convergence.
No SM electroweak completion, 1<->2 Bose, FRW or 86-eV relic.
"""
from __future__ import annotations
import math
import numpy as np
from fermion_finite_bath_pairs import _fermi_entropy
from paired_energy_quartets import implicit_one_quartet


def network_invariants(g,fc,fe):
    E=g["Echi"];k=g["Ee"];Vc=g["Vchi"];Ve=g["Ve"]
    return {
      "Nchi":float(np.dot(Vc,fc)),"Ne":float(np.dot(Ve,fe)),
      "Echi":float(np.dot(Vc*fc,E)),"Ee":float(np.dot(Ve*fe,k)),
      "S":float(sum(v*_fermi_entropy(float(f)) for v,f in zip(Vc,fc))
             +sum(v*_fermi_entropy(float(f)) for v,f in zip(Ve,fe)))
    }


def edges_from_paired_W(g):
    """Unique energy-nontrivial reversible quartets i,a ↔ j,b with i<j."""
    W=g["W"];nc,ne=len(g["Echi"]),len(g["Ee"])
    edges=[]
    for i in range(nc):
        for j in range(i+1,nc):
            for a in range(j-i,ne):
                b=a+i-j
                if W[i,a,j,b]>0:
                    edges.append((i,a,j,b))
    return edges


def evolve(g,fc0,fe0,dt,substeps=1,reverse=False):
    fc=np.asarray(fc0,dtype=float).copy()
    fe=np.asarray(fe0,dtype=float).copy()
    if (fc.shape!=g["Echi"].shape or fe.shape!=g["Ee"].shape
        or np.any(~np.isfinite(fc)) or np.any(~np.isfinite(fe))
        or np.any((fc<0)|(fc>1)) or np.any((fe<0)|(fe>1))
        or type(substeps) is not int or substeps<1
        or not math.isfinite(dt) or dt<0):
        raise ValueError("valid Pauli occupations and nonnegative time steps")
    before=network_invariants(g,fc,fe)
    edges=edges_from_paired_W(g)
    if reverse:edges=edges[::-1]
    total_abs_transfer=0.
    for _ in range(substeps):
        for i,a,j,b in edges:
            r=implicit_one_quartet(g,i,a,j,b,float(fc[i]),float(fe[a]),
                                   float(fc[j]),float(fe[b]),dt/substeps)
            fc[i],fe[a],fc[j],fe[b]=r["fi"],r["ga"],r["fj"],r["gb"]
            total_abs_transfer+=abs(r["transfer_number_GeV3"])
            if (min(fc[i],fc[j],fe[a],fe[b])<-1e-12 or
                max(fc[i],fc[j],fe[a],fe[b])>1+1e-12):
                raise ArithmeticError("discrete Fermi Pauli bound failed")
    after=network_invariants(g,fc,fe)
    return {
       "fchi":fc,"fe":fe,"edges":len(edges),"substeps":substeps,
       "Nchi_residual":after["Nchi"]-before["Nchi"],
       "Ne_residual":after["Ne"]-before["Ne"],
       "total_energy_residual":after["Echi"]+after["Ee"]
           -before["Echi"]-before["Ee"],
       "delta_Echi":after["Echi"]-before["Echi"],
       "delta_Ee":after["Ee"]-before["Ee"],
       "entropy_gain":after["S"]-before["S"],
       "abs_pair_number_flow":total_abs_transfer,
       "is_original_unmodified_5D_kernel":False
    }
