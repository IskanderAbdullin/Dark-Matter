"""Stoichiometric invariant audit for on-shell Zprime->nu+antinu event lattices.

Each event column is (-1 vector_j, +1 nu_i, +1 antinu_l).
Physical left-nullspace must include C=Nnu+Nantinu+2 NV,
L=Nnu-Nantinu, and E=∑ E_j Nj. Extra null directions reveal
ARTIFICIALLY frozen populations and possibly false equilibria.
Rank checked by independent float SVD and exact mod-prime elimination;
no sympy dependency in repository CI.
"""
from __future__ import annotations
import json
import numpy as np
from kinetics.vector_neutrino_reversible_event_graph import graph as integer_graph
from kinetics.vector_neutrino_shell_cell_overlap import graph_overlap


def stoichiometry(g):
    nv=len(g["E"])
    nk=len(g["k"])
    S=np.zeros((nv+2*nk,len(g["edges"])),dtype=np.int64)
    for col,(v,i,l,weight) in enumerate(g["edges"]):
        if weight<=0:
            raise ValueError("reaction weights must remain positive")
        S[v,col]=-1
        S[nv+i,col]=1
        S[nv+nk+l,col]=1
    return S


def exact_mod_prime_rank(matrix,prime=1000003):
    """Finite-field Gaussian elimination, independent from floating SVD."""
    A=np.array(matrix,dtype=np.int64,copy=True)%prime
    m,n=A.shape
    rank=0
    for col in range(n):
        idx=np.flatnonzero(A[rank:,col])
        if idx.size==0:continue
        pivot=rank+int(idx[0])
        if pivot!=rank:
            A[[rank,pivot],:]=A[[pivot,rank],:]
        inverse=pow(int(A[rank,col]),-1,prime)
        A[rank,col:]=(A[rank,col:]*inverse)%prime
        rows=np.flatnonzero(A[rank+1:,col])+rank+1
        if len(rows):
            factors=A[rows,col].copy()
            A[rows,col:]=(A[rows,col:]-factors[:,None]*A[rank,col:])%prime
        rank+=1
        if rank==m:break
    return rank


def audit(m=.4,h=.2,emax=4.,*,shifted=True):
    g=(graph_overlap(m,h,emax) if shifted else
       integer_graph(m,h,emax))
    S=stoichiometry(g)
    nv=len(g["E"]);nk=len(g["k"])
    C=np.concatenate((np.full(nv,2.),np.ones(nk),np.ones(nk)))
    L=np.concatenate((np.zeros(nv),np.ones(nk),-np.ones(nk)))
    E=np.concatenate((g["E"],g["k"],g["k"]))
    residuals={name:float(np.max(abs(S.T@vec)))
               for name,vec in (("C",C),("L",L),("E",E))}
    rank1=exact_mod_prime_rank(S,1000003)
    rank2=exact_mod_prime_rank(S,1000033)
    singular=np.linalg.svd(S.astype(float),compute_uv=False)
    floatrank=int(np.sum(singular>1e-9))
    isolated=np.flatnonzero(~np.any(S!=0,axis=1))
    return {"mass":m,"h":h,"emax":emax,"shifted":shifted,
            "vector_modes":nv,"nu_modes_each":nk,
            "events":S.shape[1],"states":S.shape[0],
            "rank_mod_1000003":rank1,"rank_mod_1000033":rank2,
            "rank_float_svd":floatrank,
            "nullity":S.shape[0]-rank1,
            "physical_invariant_residuals":residuals,
            "isolated_nu_indices":[int(i-nv) for i in isolated
                                   if nv<=i<nv+nk],
            "isolated_antinu_indices":[int(i-nv-nk) for i in isolated
                                       if nv+nk<=i],
            "extra_exact_invariants_beyond_C_L_E":
               S.shape[0]-rank1-3,
            "no_actual_cosmological_equilibrium_claim":True}


if __name__=="__main__":
    rows=[]
    for h in (.4,.2,.1):
        rows.append(audit(h=h,shifted=False))
        rows.append(audit(h=h,shifted=True))
    print(json.dumps(rows,indent=2))
