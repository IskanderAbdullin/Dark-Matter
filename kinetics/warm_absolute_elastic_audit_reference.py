"""Check absolute versus centered elastic rates and a coordinate-overflow stress."""
import argparse
import json
from pathlib import Path

import numpy as np

from compact_elastic_panel_reference import compact_elastic
from joint_resonant_frw_split_reference import pack
from warm_comoving_parent_competing_reference import absolute_elastic,StationaryParentElastic
from warm_comoving_parent_chart_reference import ComovingCollision,common_rule
from warm_light_occupation_stage_audit_reference import digest,read


def run(path):
    r=read(path)
    s=r["snapshot"]
    assert digest(s)==r["hash"]
    p=np.asarray(s["lab_momentum_edges_GeV"])
    Q,w=common_rule(p,(1.,),2)
    k=ComovingCollision(p,1.,Q,w)
    T=s["T_visible_GeV"]
    H=2.6e-23
    op=StationaryParentElastic(k,H,T)
    u=pack([s["entropy_variable_"+name] for name in ("scalar","fermion")])
    ref,_,t=op.reference(u,T)
    assert t==op.t
    du=.12*np.sin(.39*np.arange(len(u))+.41)
    C,dC=compact_elastic(u-ref,t,op.eta,H,op.elastic,op.d,True)
    A,dA=absolute_elastic(u,H,op.elastic,op.d)
    j0,j1=dC(du),dA(du)
    D=A.reshape(2,op.n,2)
    terms=[(c.mass+c.low)*D[j,:,0]+c.dk*D[j,:,1] for j,c in enumerate(k.c)]
    counts=np.max(abs(D[:,:,0].sum(axis=1)))/np.sum(abs(D[:,:,0]))
    energy=abs(sum(x.sum() for x in terms))/sum(abs(x).sum() for x in terms)
    assert max(counts,energy)<3e-12
    def rel(x,y):
        return float(np.linalg.norm(x-y)/np.linalg.norm(y))
    eps=2e-5
    df=(absolute_elastic(u+eps*du,H,op.elastic,op.d)[0]-
        absolute_elastic(u-eps*du,H,op.elastic,op.d)[0])/(2*eps)
    metrics=[rel(C,A),rel(j0,j1),rel(df,j1)]
    assert metrics[0]<4e-11 and metrics[1]<4e-11 and metrics[2]<5e-7
    # Same ordinary finite occupation, represented relative to an
    # deliberately cold synthetic reference. Only the coordinate changes.
    t0=1e-5
    artificial=pack([np.array([(c.mass+c.low)/t0,c.dk/t0]) for c in k.c])
    caught=False
    with np.errstate(over="raise", invalid="raise"):
        try:
            compact_elastic(u-artificial,t0,(0.,0.),H,op.elastic,op.d,True)
        except FloatingPointError:
            caught=True
    assert caught
    return dict(input_hash=r["hash"],panels=len(op.elastic),centered_vs_absolute_flux=metrics[0],
                centered_vs_absolute_tangent=metrics[1],absolute_vs_finite_difference=metrics[2],
                count_and_energy_relative_residuals=[float(counts),float(energy)],
                same_physical_state_under_artificial_cold_reference_requires_absolute=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("state",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.state),indent=2,allow_nan=False))
