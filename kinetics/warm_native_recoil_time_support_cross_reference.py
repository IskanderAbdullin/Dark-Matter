"""Compare two native time ladders before the corrected cold-shell continuation."""

import argparse
import json
from pathlib import Path

import numpy as np

from frw_exterior_restart_reference import import_restart
from joint_prompt_two_moment_evolution_reference import Collision
from late_radiation_background_reference import MPL,TabulatedRadiation
from three_species_decay_onsager_reference import mp
from warm_coupled_serial_restart_reference import sha
from warm_history_lower_matched_audit_reference import measures
from warm_light_occupation_stage_audit_reference import digest,read
from warm_near_rest_light_inverse_decay_reference import window


def run(eos,old,fine):
    def read_at(D,i):
        q=read(D/f"step-{i}.json")
        assert digest(q["snapshot"])==q["snapshot_sha256"]
        p=q["plan"]
        assert (p.get("EOS_sha256") or p.get("eos_sha256"))==sha(eos)
        return q
    for q in (32,64):
        a=read_at(old/f"q{q}",8)
        s=read_at(fine/f"q{q}",1)
        assert (s.get("previous") or s.get("parent_snapshot_sha256"))==a["snapshot_sha256"]
    A=read_at(old/"q32",8)["snapshot"]
    p0=np.asarray(A["lab_momentum_edges_GeV"])
    p0=p0[:-1]
    n=len(p0)
    scale=measures(A).sum(axis=-1)*A["scale_factor"]**3
    result=[]
    for i in range(9):
        paths=[(old/"q32",i+8),(old/"q64",i+8),
               (fine/"q32",2*i),(fine/"q64",2*i)]
        if i==0: paths[2:]=paths[:2]
        rows=[read_at(D,j) for D,j in paths]
        T=rows[0]["snapshot"]["T_visible_GeV"]
        R=TabulatedRadiation(eos,T)
        hist=[]; moment=[]
        edge=None
        for r in rows:
            S=r["snapshot"]
            np.testing.assert_allclose(T,S["T_visible_GeV"],rtol=1.5e-13,atol=0.)
            p=np.asarray(S["lab_momentum_edges_GeV"])
            if edge is None: edge=p[:n+1]
            np.testing.assert_allclose(edge,p[:n+1],rtol=1.5e-13,atol=0.)
            W=measures(S)
            E=W[1].sum()+import_restart(S["positive_exterior_restart"]).energy(1.)
            H=np.sqrt((R.state(1.)[1]+E)/3)/MPL
            op=Collision(p)
            th=[np.asarray(S["entropy_variable_"+k]) for k in ("scalar","fermion")]
            active=int(np.sum(op.J==0))
            V=op.c[1].vol[0]
            rate=float(op.rate(th)[1,0,0]/V/H)
            closed=bool(np.hypot(mp,p[-1])<window(p[1])[-2])
            if closed: assert active==0 and rate==0.
            hist.append(dict(snapshot_sha256=r["snapshot_sha256"],first_occupation=float(W[0,1,0]/V),
                             geometry_closed=closed,active_prompt_nodes=active,first_prompt_rate_over_H=rate))
            moment.append(W[:,:,:n]*S["scale_factor"]**3)
        a,b,c,d=moment
        def L(x):
            return (np.sum(abs(x),axis=-1)/scale).tolist()
        row=dict(coarse_frame=i,T_GeV=T,relative_a=rows[0]["snapshot"]["scale_factor"]/A["scale_factor"],
                 branches=hist,time_change_q32=L(c-a),time_change_q64_shared=L(d-b),
                 support_change_coarse=L(b-a),support_change_fine=L(d-c),interaction=L((d-b)-(c-a)))
        if i:
            dt=np.log(row["relative_a"]/result[-1]["relative_a"])
            row["cold_slopes_shared_intervals"]=[(x["first_occupation"]-y["first_occupation"])/dt
                                                 for x,y in zip(hist,result[-1]["branches"],strict=True)]
        result.append(row)
    return dict(scope="Only the native before-correction recoil interval and matched shared cells. "
                      "The plotted geometric gap is unaffected by time interpolation of spectra.",
                metric_order=["number","energy","pressure"],branch_order=["coarse q32","coarse q64","fine q32","fine q64"],frames=result)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","coarse","fine"):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.coarse,a.fine),indent=2,allow_nan=False))
