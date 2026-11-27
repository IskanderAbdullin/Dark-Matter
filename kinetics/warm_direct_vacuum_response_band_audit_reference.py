"""Distinguish the rare equality-soft tail from the bulk scalar tracer response."""
import argparse,json
from pathlib import Path
import numpy as np
from archived_collisionless_streaming_reference import Background
from archived_scalar_pair_reaction_reference import points
from finite_cooling_relativistic_los_reference import trajectory
from joint_visible_heavy_adiabatic_reference import snapshot_functions
from three_species_decay_onsager_reference import mp
from warm_direct_vacuum_fixed_metric_response_reference import kernel,layout
from warm_light_occupation_stage_audit_reference import digest,read


def run(eos,archive,source,metrics):
    bg=Background(eos)
    A=read(archive)
    S=A["snapshot"]
    assert digest(S)==A["snapshot_sha256"]
    c=read(source)
    e=np.array([x if x is not None else 36. for x in c["conserved_momentum_edges_q_over_T0"]])
    m=np.asarray(c["controls"][-1]["shell_yields"]).sum(axis=0)
    j=int(np.flatnonzero(e==.1)[0])
    G=m/m.sum()
    T=S["T_visible_GeV"]
    a=bg.scale(T)
    f,*_=snapshot_functions(S)
    p=np.asarray(S["lab_momentum_edges_GeV"])
    p=np.unique(np.r_[p,.1*bg.T0/a])
    rows=[]
    for file in metrics:
        mt=trajectory(file)
        s=[]
        for n,time in ((8,7500),(12,13000)):
            q,w=layout(e,n,"q")
            p0,_,h=points(p,mp,f,n)
            qh=a*p0/bg.T0
            eh,_,_,R=kernel(bg,mt,np.r_[q,qh],time)
            sE=eh[:len(q)].reshape(len(G),n)@w
            sU=(eh[:len(q)]*R[:len(q)]).reshape(len(G),n)@w
            he=eh[len(q):]*h
            hu=he*R[len(q):]
            col=qh<.1
            coarse=np.array([[float((G*sU).sum()),float((G*sE).sum()),float((G[:j]*sU[:j]).sum()),float((G[:j]*sE[:j]).sum())],
                             [float(hu.sum()),float(he.sum()),float(hu[col].sum()),float(he[col].sum())]])
            whole=coarse[:,0]/coarse[:,1]
            soft=coarse[:,2]/coarse[:,3]
            remainder=(coarse[:,0]-coarse[:,2])/(coarse[:,1]-coarse[:,3])
            s.append(dict(q_order=n,time_nodes=time,whole_response=whole.tolist(),cold_only_response=soft.tolist(),
                         cold_energy_weight_fraction=(coarse[:,3]/coarse[:,1]).tolist(),
                         not_cold_only_response=remainder.tolist(),
                         direct_minus_held_whole=float(whole[0]-whole[1]),
                         direct_minus_held_not_cold=float(remainder[0]-remainder[1]),
                         cold_removal_change_in_shape_contrast=float((remainder[0]-remainder[1])-(whole[0]-whole[1]))))
        assert abs(s[-1]["cold_removal_change_in_shape_contrast"])<.001
        rows.append(dict(k_per_Mpc=mt[0],controls=s))
    return dict(scope="Linear functional decomposition into q/T0≤.1 and its complement on fixed driver, "
                      "not an independently generated metric or solved production perturbation.",
                shape_order=["direct vacuum mix","held q64"],momentum_cut_over_T0=.1,rows=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","archive","source"):
        p.add_argument(k,type=Path)
    p.add_argument("metrics",nargs="+",type=Path)
    a=p.parse_args();print(json.dumps(run(a.eos,a.archive,a.source,a.metrics),indent=2,allow_nan=False))
