"""Direct-occupation prompt solve for a high-contrast known-cohort stress test."""

import argparse
import json
from pathlib import Path

import numpy as np

from joint_prompt_two_moment_evolution_reference import Collision
from joint_resonant_frw_split_reference import pack
from late_archived_shape_continuation_reference import Late,collide
from late_radiation_background_reference import MPL,TabulatedRadiation
from warm_cohort_two_packet_geometry_reference import highrate,mom,seed
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import pad


class DirectPrompt(Late):
    def moment(self,q):
        scalar=self.arrays(q)[0]
        # Same physical domain, with a numerical guard far outside all
        # admitted stress occupations so the analytic f^2 metric stays finite.
        if np.min(scalar[0]+np.minimum(0.,scalar[1]))<=1e-130:
            raise ValueError("invalid scalar endpoint for direct solve")
        return super().moment(q)

    def centered(self,z,ref,R0):
        assert self.prompt and not self.panels
        return self.collision(ref+z)


def tests(S,H):
    q=np.asarray(S["lab_momentum_edges_GeV"])
    th=pack(pad(S,q))
    grid=Collision(q)
    coarse=Late(grid,H,[],True,local_jac=True)
    direct=DirectPrompt(grid,H,[],True)
    ref,R0,_=coarse.reference(th,S["T_visible_GeV"])
    f,j=coarse.centered(th-ref,ref,R0)
    g,k=direct.centered(th-ref,ref,R0)
    v=.2*np.sin(.51*np.arange(len(th))+.3)
    x,y=j(v),k(v)
    e=2e-6
    ytest=(direct.centered(th+e*v-ref,ref,R0)[0]-direct.centered(th-e*v-ref,ref,R0)[0])/(2*e)
    errors=[np.linalg.norm(f-g)/np.linalg.norm(f),np.linalg.norm(x-y)/np.linalg.norm(x),
            np.linalg.norm(y-ytest)/np.linalg.norm(y)]
    assert max(errors[:2])<2e-9 and errors[-1]<2e-7
    return [float(z) for z in errors]


def evolve(S,hard,eos,split):
    A,ledger=seed(S,hard)
    p=np.asarray(A["lab_momentum_edges_GeV"])
    keys=("entropy_variable_scalar","entropy_variable_fermion")
    U=pack([A[k] for k in keys])
    n,k,E=mom(A)
    H=np.sqrt((TabulatedRadiation(eos,S["T_visible_GeV"]).state(1.)[1]+E.sum())/3)/MPL
    out=[]
    init=highrate(A,1.,eos)
    base=n.sum(axis=1)
    M=DirectPrompt(Collision(p),H,[],True)
    times=[0.,1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5]
    for a,b in zip(times[:-1],times[1:],strict=True):
        r=0.; accepted=rejected=0
        for _ in range(split):
            U,v,acc,rej=collide(M,U,(b-a)/split,S["T_visible_GeV"],continuation=True,balanced=True)
            r=max(r,v); accepted+=acc; rejected+=rej
        st=dict(A,**{k:z.tolist() for k,z in zip(keys,M.arrays(U),strict=True)})
        N,K,E1=mom(st)
        errors=[float(np.max(abs(N.sum(axis=1)/base-1.))),float(abs(E1.sum()/E.sum()-1.))]
        assert max(errors)<2e-12
        out.append(dict(H_times_elapsed=b,high_inverse_and_tagged_energy_rate=highrate(st,1.,eos),
                        solver=[r,accepted,rejected],relative_number_and_energy_errors=errors))
        print(json.dumps(dict(hard=hard,subdivision=split,t=b)),file=__import__("sys").stderr,flush=True)
    return dict(hard=hard,seed_ledger=ledger,initial=init,controls=out)


def run(eos,D):
    A=read(D/"warm-light-occupation-next-q32-tape-2026-11-26.json.gz")[0]
    S=A["snapshot"]
    assert A["snapshot_sha256"]==digest(S)
    E=mom(S)[2].sum()
    H=np.sqrt((TabulatedRadiation(eos,S["T_visible_GeV"]).state(1.)[1]+E)/3)/MPL
    check=tests(S,H)
    rows=[evolve(S,p,eos,n) for p in (.52,.58) for n in (1,2)]
    return dict(scope="Controlled fixed-clock prompt-only stress. Direct occupation formula avoids "
                      "forming the prior exponentially large trial factors. Not a general solver guarantee.",
                ordinary_archive_flux_JVP_finite_difference_defects=check,cases=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path);p.add_argument("--tapes",type=Path,default=Path("results"))
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.tapes),indent=2,allow_nan=False))
