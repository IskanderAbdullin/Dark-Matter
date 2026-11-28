"""Short number- and energy-closed prompt transient with reacting exterior."""

import argparse
import json
from pathlib import Path

import numpy as np

from joint_prompt_two_moment_evolution_reference import Collision
from joint_resonant_frw_split_reference import pack
from late_archived_shape_continuation_reference import Late,collide
from three_species_decay_onsager_reference import width
from warm_light_occupation_stage_audit_reference import digest,read


def transient(source,grid,ncore,H,subdivision):
    th=[np.asarray(source["entropy_variable_"+s]) for s in ("scalar","fermion")]
    both=np.array([np.pad(t,((0,0),(0,len(grid)-ncore-1))) for t in th])
    both[:,0,ncore:]=1000.  # zero to double precision in both statistics
    model=Late(Collision(grid),H,[],True,local_jac=True)
    u=pack(both)
    def evaluation(u):
        C=model.moment(u)[0].reshape(2,model.n,2)
        N=C[:,:,0].sum(axis=1)
        E=sum(np.sum(C[j,:,0]*(d.mass+d.low)+C[j,:,1]*d.dk) for j,d in enumerate(model.k.c))
        return N,E,C[:,ncore:,0].sum(axis=1),C[:,-1,0]
    A=evaluation(u)
    assert np.all(A[2]==0.)
    times=(0.,1e-8,5e-8,2e-7,1e-6,5e-6)
    rows=[]
    for p,q in zip(times[:-1],times[1:],strict=True):
        residual=0.;acc=rej=0
        for _ in range(subdivision):
            u,r,a,b=collide(model,u,(q-p)/subdivision,source["T_visible_GeV"],
                            continuation=True,balanced=True)
            residual=max(residual,r);acc+=a;rej+=b
        N,E,X,L=evaluation(u)
        errors=[float(np.max(abs(N/A[0]-1))),float(abs(E/A[1]-1))]
        assert max(errors)<1e-10
        rows.append(dict(elapsed_times_H=q,central_solve_residual=residual,
                         accepted=acc,rejected=rej,relative_number_and_energy_defects=errors,
                         beyond_original_edge_fraction=(X/A[0]).tolist(),
                         outermost_shell_fraction=(L/A[0]).tolist()))
    return rows


def run(D):
    fits=[read(D/f"warm-light-occupation-next-q{q}-tape-2026-11-26.json.gz")[0] for q in (16,32)]
    A,B=[v["snapshot"] for v in fits]
    assert all(v["snapshot_sha256"]==digest(S) for v,S in zip(fits,(A,B),strict=True))
    P,Q=[np.asarray(s["lab_momentum_edges_GeV"]) for s in (A,B)]
    assert np.allclose(Q[:len(P)],P,rtol=3e-14,atol=0.)
    H=read(D/"warm-near-rest-light-inverse-decay-2026-11-26.json.gz")["cases"][0]["stages"][0]["fine"]["H_over_rest_decay_width"]*width
    closed=transient(A,P,len(P)-1,H,1)
    assert np.all(np.asarray(closed[-1]["beyond_original_edge_fraction"])==0.)
    a=transient(A,Q,len(P)-1,H,1)
    b=transient(A,np.r_[Q,Q[-1]*np.geomspace(1.,2.,6)[1:]],len(P)-1,H,1)
    refined=transient(A,Q,len(P)-1,H,2)
    x,y,z=[np.asarray(t[-1]["beyond_original_edge_fraction"]) for t in (a,b,refined)]
    xerr=float(np.max(abs(x/y-1.))); terr=float(np.max(abs(x/z-1.)))
    assert max(xerr,terr)<.006
    ref=read(D/"warm-post-checkpoint-open-prompt-flux-2026-11-26.json.gz")["cases"][0]["stages"][0]
    linear=a[-1]["elapsed_times_H"]*np.asarray(ref["full_support"][-1]["rate_over_current_species_count_per_Hubble"])
    return dict(scope="Frozen radiation clock and energy background, explicit positive vacuum exteriors "
                      "with reversible represented prompt transfer only. No expansion, new sources, elastic "
                      "collision, physical inherited exterior, or late evolution is inferred.",
                q16_checkpoint_sha256=fits[0]["snapshot_sha256"],q16_original=closed,q32_padded=a,
                q64_padded=b,q32_halved_steps=refined,
                relative_terminal_support_and_time_changes=[xerr,terr],
                frozen_linear_no_reentry_fractions_at_same_time=linear.tolist(),
                ratio_of_frozen_linear_to_q32_population_end=(linear/x).tolist())


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--tapes",type=Path,default=Path("results"))
    print(json.dumps(run(p.parse_args().tapes),indent=2,allow_nan=False))
