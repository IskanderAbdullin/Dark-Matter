"""Detect non-physical collision invariants before trusting adaptive kinetic grids.

For every active resonant reaction chi1_k <-> phi_i + chi2_j, a
collision invariant (b_i,c_j,m_k) must obey b_i+c_j=m_k.

Eliminating each parent variable gives constraints
 b_i+c_j - b_i0-c_j0 = 0
for all pairs sharing one parent energy. Their nullspace includes
exactly three EXPECTED quantities: boson number offset, daughter-fermion
number offset, and energy (b_i=c_i=e_i). Any additional dimensions
are possible SPURIOUS conservation laws of the finite collision graph.

The explicit physical regime can also have disconnected thresholds; the
rank counts algebraic invariants, not an observational proof of mixing.
Reference: Bernhoff, Discrete Quantum Kinetic Equation,
La Matematica 2 (2023) 836-860, doi:10.1007/s44007-023-00070-1.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from scipy.sparse import coo_matrix


def three_block_indices():
    return list(range(16))+list(range(2500,2520))+list(range(5000,5020))


def sparse_cold_hot_indices():
    return list(range(8))+[12,16,24,36,50,75,100,150,200,350,500,
                           750,1000,1500,2000,3000,4000,5000,6000]


def audit_integer_grid(indices,delta=.001,mass=.1,declared_parent_keys=None):
    index=np.asarray(sorted(set(int(i) for i in indices)),dtype=np.int64)
    if len(index)<3 or np.any(index<0) or delta<=0 or mass<=0:
        raise ValueError("Invalid resonant daughter-grid indices/parameters")
    n=len(index)
    energies=(index+.5)*delta
    pairs=np.argwhere(4*energies[:,None]*energies[None,:]>=mass*mass)
    if not len(pairs):
        raise ValueError("No physical 1<->2 pair on this grid")
    parent_keys=index[pairs[:,0]]+index[pairs[:,1]]+1
    active_parent=np.unique(parent_keys)
    grouped={}
    for (i,j),parent_key in zip(pairs,parent_keys):
        grouped.setdefault(int(parent_key),[]).append((int(i),int(j)))
    rows=[];cols=[];values=[];count=0
    for pp in grouped.values():
        i0,j0=pp[0]
        for i,j in pp[1:]:
            for col,weight in ((i,1.),(n+j,1.),(i0,-1.),(n+j0,-1.)):
                rows.append(count);cols.append(col);values.append(weight)
            count+=1
    constraints=coo_matrix((values,(rows,cols)),shape=(count,2*n)).tocsr()
    gram=(constraints.T@constraints).toarray()
    eigen,vec=eigh(gram)
    threshold=max(1e-9,float(max(eigen))*1e-11)
    rank=int(np.count_nonzero(eigen>threshold))
    nullity=int(2*n-rank)

    # Three analytic invariants of the continuum on-shell collision
    # process. Verify that all lie in this discrete nullspace.
    phys=np.zeros((2*n,3))
    phys[:n,0]=1.
    phys[n:,1]=1.
    phys[:n,2]=energies
    phys[n:,2]=energies
    analytic_residual=float(np.max(abs(constraints@phys))) if count else 0.
    q,_=np.linalg.qr(phys,mode="reduced")
    witness=None
    if nullity>3:
        candidates=vec[:,eigen<=threshold]
        projected=candidates-q@(q.T@candidates)
        norms=np.linalg.norm(projected,axis=0)
        k=int(np.argmax(norms))
        if norms[k]>1e-8:
            w=projected[:,k]/norms[k]
            witness={
                "extra_invariant_max_reaction_constraint_residual":float(
                    np.max(np.abs(constraints@w))) if count else 0.,
                "orthogonality_to_three_physical_invariants":float(
                    np.max(np.abs(q.T@w))),
                "max_absolute_coefficient":float(np.max(np.abs(w))),
            }
    inactive_parent=[]
    if declared_parent_keys is not None:
        declared=set(int(x) for x in declared_parent_keys)
        inactive_parent=sorted(declared-set(map(int,active_parent)))
        if set(active_parent)-declared:
            raise ValueError("Active parent not in declared grid")
    active_b=set(map(int,pairs[:,0]))
    active_f=set(map(int,pairs[:,1]))
    return {
        "daughter_grid_n":n,
        "daughter_delta":float(delta),
        "daughter_energy_min":float(energies[0]),
        "daughter_energy_max":float(energies[-1]),
        "mother_mass_over_T":float(mass),
        "active_pair_reactions":int(len(pairs)),
        "active_unique_mother_energy_modes":int(len(active_parent)),
        "linear_constraints_after_parent_elimination":int(count),
        "rank_on_daughter_coefficients":rank,
        "active_graph_collision_invariants":nullity,
        "expected_physical_invariants":3,
        "additional_invariants_beyond_number_number_energy":nullity-3,
        "analytical_three_invariant_residual":analytic_residual,
        "inactive_declared_parent_modes":len(inactive_parent),
        "inactive_declared_parent_keys":inactive_parent,
        "full_grid_collision_invariants_including_inactive_parents":
            nullity+len(inactive_parent),
        "isolated_boson_modes":n-len(active_b),
        "isolated_chi2_modes":n-len(active_f),
        "extra_invariant_witness":witness,
        "accept_for_unique_chemical_equilibrium":
            nullity==3 and not inactive_parent and len(active_b)==n and len(active_f)==n,
        "caveat":"Rank is algebraic connectivity of sampled reactions, not a claim that continuum phase-space kinetics mixes all modes. Additional invariants can distort equilibria on sparse grids.",
    }



def extra_invariant_coefficients(indices,delta=.001,mass=.1):
    """Construct ONE nonphysical left collision nullvector for an adaptive graph.

    The returned parent keys match the unique active sums used by
    AdaptiveClosedDecay. This is a DIRECT invariance witness for the
    full nonlinear ODE, not only its equilibrium linearization.
    """
    idx=np.asarray(sorted(set(int(i) for i in indices)),dtype=np.int64)
    n=len(idx);energy=(idx+.5)*delta
    i,j=np.nonzero(4*energy[:,None]*energy[None,:]>=mass*mass)
    keys=idx[i]+idx[j]+1
    unique=np.unique(keys)
    grouped={}
    for a,b,k in zip(i,j,keys):
        grouped.setdefault(int(k),[]).append((int(a),int(b)))
    rows=[];cols=[];vals=[];z=0
    for pairs in grouped.values():
        ai,aj=pairs[0]
        for bi,bj in pairs[1:]:
            for col,v in ((bi,1.),(n+bj,1.),(ai,-1.),(n+aj,-1.)):
                rows.append(z);cols.append(col);vals.append(v)
            z+=1
    C=coo_matrix((vals,(rows,cols)),shape=(z,2*n)).tocsr()
    eig,V=eigh((C.T@C).toarray())
    tol=max(1e-9,float(max(eig))*1e-11)
    phys=np.zeros((2*n,3))
    phys[:n,0]=1;phys[n:,1]=1
    phys[:n,2]=energy;phys[n:,2]=energy
    Q,_=np.linalg.qr(phys,mode="reduced")
    projected=V[:,eig<=tol]-Q@(Q.T@V[:,eig<=tol])
    norms=np.linalg.norm(projected,axis=0)
    best=int(np.argmax(norms))
    if norms[best]<1e-8:
        raise ValueError("No nonphysical invariant found on this graph")
    w=projected[:,best]/norms[best]
    parent=np.array([w[grouped[int(k)][0][0]]+
                     w[n+grouped[int(k)][0][1]] for k in unique])
    residue=max(abs(C@w)) if z else 0.
    assert residue<1e-8 and np.max(abs(Q.T@w))<1e-8
    return {"boson_coefficients":w[:n],
            "chi2_coefficients":w[n:],
            "parent_coefficients":parent,
            "parent_integer_energy_keys":unique,
            "max_reaction_residual":float(residue)}


def run():
    cases={}
    for nq,delta in ((30,.08),(50,.08),(30,.02),(30,.005)):
        declared=np.arange(1,2*nq,dtype=int)
        cases[f"uniform_n{nq}_delta{delta}"]=audit_integer_grid(
            range(nq),delta=delta,declared_parent_keys=declared)
    cases["sparse_cold_hot_27"]=audit_integer_grid(sparse_cold_hot_indices())
    cases["three_blocks_56"]=audit_integer_grid(three_block_indices())
    for k,v in cases.items():
        assert v["analytical_three_invariant_residual"]<1e-10,(k,v)
        if v["extra_invariant_witness"] is not None:
            assert v["extra_invariant_witness"]["extra_invariant_max_reaction_constraint_residual"]<1e-8,(k,v)
    return {
        "scope":"Independent algebraic rank analysis, no ODE, no CLASS or observations",
        "reference":"Bernhoff 2023 doi:10.1007/s44007-023-00070-1; finite reaction-graph approach independent of author's numerical model",
        "physical_count":3,
        "cases":cases,
        "warning":"Uniform historic grids have extra conserved inactive mother slots; sparse adaptive candidate graphs can contain additional nonphysical invariants even after removing inactive mother slots."
    }


if __name__=="__main__":
    pa=argparse.ArgumentParser()
    pa.add_argument("--output",default="outputs/collision_invariant_rank.json")
    args=pa.parse_args()
    data=run()
    dst=Path(args.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(data,indent=2),encoding="utf-8")
    print(json.dumps(data,indent=2))
