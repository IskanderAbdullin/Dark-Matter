"""Measure the homogeneous stable-fermion/CDM difference in the matched reference."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from archived_collisionless_streaming_reference import Background
from archived_scalar_pair_reaction_reference import points
from joint_visible_heavy_adiabatic_reference import snapshot_functions
from three_species_decay_onsager_reference import m
from warm_light_occupation_stage_audit_reference import digest,read


def run(eos,archive):
    h=.6736
    b=Background(eos,omega_m=(.12+.02237)/h**2)
    ref=read(archive);S=ref["snapshot"]
    assert ref["snapshot_sha256"]==digest(S)
    T=S["T_visible_GeV"]
    a0=b.scale(T)
    e=np.asarray(S["lab_momentum_edges_GeV"])
    _,_,F,_=snapshot_functions(S)
    omega=.12/h**2
    dE=lambda x:x*x/(np.hypot(x,1.)+1.)
    stamps=np.array([a0,1e-10,1e-8,1e-6,b.aeq,1/21,1.])
    reports=[]
    for n in (16,28):
        p,_,w=points(e,m,F,n);q=p*a0/m;w/=sum(w)
        tod=float(w@dE(q))
        def dh(a):
            add=(dE(q[None,:]/a[:,None])@w-tod)/(1+tod)
            add[a==1.]=0.
            E=b.weighted_radiation(a)+b.om*a+b.ol*a**4
            ratio=omega*a*add/E
            return ratio/(np.sqrt(1+ratio)+1)
        l=np.unique(np.r_[np.log(a0),np.log(b.aeq),np.clip(b.v,np.log(a0),np.log(b.aeq)),
                           np.linspace(np.log(a0),np.log(b.aeq),30)])
        z,v=roots_legendre(8)
        s=np.diff(l)[:,None]/2
        a=np.exp((l[:-1,None]+l[1:,None])/2+s*z).ravel()
        wtime=(s*v).ravel()*a
        E=np.sqrt(b.weighted_radiation(a)+b.om*a+b.ol*a**4)
        effect=float(-b.D0*((wtime/E)@(dh(a)/(1+dh(a)))))
        reports.append(dict(nodes_per_cell=n,homogeneous_delta_H_over_H=dh(stamps).tolist(),
                            delta_conformal_Mpc_from_anchor_to_equality=effect))
    np.testing.assert_allclose(reports[0]["homogeneous_delta_H_over_H"],reports[1]["homogeneous_delta_H_over_H"],rtol=2e-9,atol=1e-27)
    return dict(scope="Only archived stable phase free redshifting after the 7.394 MeV anchor, rescaled to the same "
                      "present dark density as CDM. Homogeneous proxy, no self-consistent production, clustering "
                      "backreaction, or total bound for unrepresented phases.",
                archive=ref["snapshot_sha256"],scales=stamps.tolist(),controls=reports)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("eos",type=Path);p.add_argument("archive",type=Path)
    a=p.parse_args();print(json.dumps(run(a.eos,a.archive),indent=2,allow_nan=False))
