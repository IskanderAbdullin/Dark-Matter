"""Bound conditional path means and quantiles without choosing a shape in each shell."""
import argparse,json
from pathlib import Path
import numpy as np
from archived_collisionless_streaming_reference import Background
from archived_scalar_pair_reaction_reference import points
from joint_visible_heavy_adiabatic_reference import snapshot_functions
from three_species_decay_onsager_reference import mp
from warm_light_occupation_stage_audit_reference import digest,read


def run(eos,archive,psd):
    b=Background(eos)
    d=read(archive);S=d["snapshot"]
    assert d["snapshot_sha256"]==digest(S)
    a=b.scale(S["T_visible_GeV"])
    P=read(psd)
    q=np.array([float(x) for x in P["conserved_momentum_edges_q_over_T0"][:-1]])
    V=np.array(P["controls"][-1]["shell_yields"])
    V=np.r_[V,V.sum(axis=0)[None,:]]
    V/=V.sum(axis=1)[:,None]
    qn=np.maximum(q,1e-30)*b.T0
    gmin=np.min(b.gr(b.u))
    gsmax=np.max(b.gs(b.u))
    rad=b.ogamma/2*np.exp(gmin+4*(np.log(b.gs0)-gsmax)/3)
    top=b.D0*(b.aeq-a)/np.sqrt(min(rad,b.or0))
    d0=b.streams(qn,mp,a,b.aeq,7);d0[0]=0.
    ref=b.streams(qn,mp,a,b.aeq,10);ref[0]=0.
    np.testing.assert_allclose(ref,d0,rtol=2e-11,atol=2e-10)
    assert np.all(np.diff(d0)>0.) and top>d0[-1]
    d1=np.r_[d0,top]
    F,*_=snapshot_functions(S)
    p,_,w=points(np.asarray(S["lab_momentum_edges_GeV"]),mp,F,24)
    X=b.streams(a*p,mp,a,b.aeq,7)
    z=w@X/w.sum()
    c=np.cumsum(V,axis=1)
    qs=np.array([.1,.5,.9])
    quant=[]
    for v in c:
        j=np.searchsorted(v,qs)
        quant.append([(float(d1[i]),float(d1[i+1])) for i in j])
    return dict(scope="Distance since the common warm archive time to equality, conditional on both sources "
                      "already deposited and collisionless. Monotone paths give shell bounds; not a transfer.",
                held_mean_path_Mpc=float(z),source_order=["lepton","measured charged pion","yield weighted two sources"],
                positive_shell_mean_distance_bounds_Mpc=np.stack([V@d1[:-1],V@d1[1:]],axis=1).tolist(),
                percentile_levels=qs.tolist(),shell_only_percentile_distance_enclosures_Mpc=quant,
                comparison_anchor_scale=float(a))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","archive","psd"):
        p.add_argument(k,type=Path)
    a=p.parse_args();print(json.dumps(run(a.eos,a.archive,a.psd),indent=2,allow_nan=False))
