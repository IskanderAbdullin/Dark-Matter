"""Verify atomic q64 guided snapshots independently of the solver wrapper."""

import argparse
import json
from pathlib import Path

import numpy as np

from frw_exterior_restart_reference import import_restart
from warm_coupled_serial_restart_reference import sha
from warm_light_occupation_stage_audit_reference import digest,read
from warm_observed_suboperator_continuation_reference import metrics
from warm_prompt_reentry_shape_audit_reference import get_mom


def run(start,guides,tape):
    A=read(start)["snapshot"]
    G=[x["checkpoint"] for x in read(guides)["blocks"]]
    states=read(tape)
    assert len(states)==len(G)
    a0=A["scale_factor"]
    def mom(s):
        return get_mom(np.asarray(s["lab_momentum_edges_GeV"]),
                       [s["entropy_variable_scalar"],s["entropy_variable_fermion"]])
    N0=mom(A)[0].sum()*a0**3
    ext0=import_restart(A["positive_exterior_restart"]).array()
    n0=ext0.shape[1]
    prev=digest(A)
    out=[]
    for i,(C,G) in enumerate(zip(states,G,strict=True),1):
        S=C["snapshot"]; X=G["snapshot"]
        assert C["step"]==i and C["parent_snapshot_sha256"]==prev
        assert C["snapshot_sha256"]==digest(S)
        assert C["plan"]["source_sha256"]==sha(start) and C["plan"]["guides_sha256"]==sha(guides)
        assert C["plan"]["guide_snapshots"][i-1]==G["snapshot_sha256"]==digest(X)
        prev=digest(S)
        now=S["scale_factor"]
        assert abs(now/X["scale_factor"]-1)<3e-14
        wide=np.asarray(S["lab_momentum_edges_GeV"]); narrow=np.asarray(X["lab_momentum_edges_GeV"])
        np.testing.assert_allclose(wide[:len(narrow)],narrow,rtol=5e-14,atol=0.)
        ext=import_restart(S["positive_exterior_restart"]).array()
        clone=ext[:,:n0].copy()
        b=now/a0
        clone[[0,5]]*=np.array([b,b**3])[:,None]
        np.testing.assert_allclose(clone,ext0,rtol=9e-14,atol=0.)
        N=mom(S)[0]
        total=now**3*(N.sum()+ext[5,n0:].sum())
        err=abs(total/N0-1.)
        assert err<2e-12
        e=C["operator_summary"]
        assert abs(e["energy_ledger"])<2e-12 and abs(e["count_ledger"])<2e-12
        m=metrics([S["entropy_variable_scalar"],S["entropy_variable_fermion"]])
        mn=metrics([X["entropy_variable_scalar"],X["entropy_variable_fermion"]])
        out.append(dict(step=i,temperature_GeV=S["T_visible_GeV"],snapshot_sha256=prev,
                        cumulative_number_defect=float(err),local_energy_work_defect=e["energy_ledger"],
                        inherited_exterior_count=n0,new_collisionless_receipts=ext.shape[1]-n0,
                        q64_archive_peak=m,q32_archive_peak=mn,
                        common_count_L1=C["audit"]["common_grid_count_L1"],
                        solver_iterations=C["audit"]["trial"].get("newton_iterations"),
                        used_guided_guess=C["audit"]["trial"]["used_q32_guess"]))
    return dict(scope="Checks accepted atomic snapshots. Reconstructed extrema do not bound any "
                      "unknown sub-shell distribution or actual heavy occupation.",steps=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for a in ("start","guides","tape"):p.add_argument(a,type=Path)
    x=p.parse_args()
    print(json.dumps(run(x.start,x.guides,x.tape),indent=2,allow_nan=False))
