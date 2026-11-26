"""Read the seven independent q64 continuation blocks and their matching q32 controls."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np

from frw_exterior_restart_reference import import_restart
from late_radiation_background_reference import MPL,TabulatedRadiation
from warm_light_occupation_stage_audit_reference import digest,read
from warm_observed_suboperator_continuation_reference import metrics
from warm_prompt_high_parent_ledger_reference import high
from warm_prompt_reentry_shape_audit_reference import get_mom


def run(eos,initial,tapes,q32=None):
    raw=initial.read_bytes()
    if initial.suffix==".gz":
        raw=gzip.decompress(raw)
    seed=json.loads(raw)
    S0=seed[-1]["snapshot"]
    prev=seed[-1]["snapshot_sha256"]
    assert digest(S0)==prev
    a0=S0["scale_factor"]
    def mom(s):
        return get_mom(np.asarray(s["lab_momentum_edges_GeV"]),
                       [s["entropy_variable_scalar"],s["entropy_variable_fermion"]])
    N0,_,E0=mom(S0)
    old=import_restart(S0["positive_exterior_restart"]).array()
    idx=old.shape[1]
    target=N0.sum()*a0**3
    frames=[]
    for path in tapes:
        frames.extend(read(path))
    assert len(frames)==7
    archive=read(q32) if q32 else None
    if archive:
        assert len(archive)==8
    rr=[]
    for i,C in enumerate(frames,1):
        assert C["index"]==i and C["previous"]==prev
        S=C["snapshot"];prev=digest(S)
        assert C["snapshot_sha256"]==prev
        assert C["plan"]["tape_sha256"]==hashlib.sha256(raw).hexdigest()
        z=S["scale_factor"]; b=z/a0
        assert abs(b/C["plan"]["step_expansion"]**i-1)<2e-14
        ext=import_restart(S["positive_exterior_restart"])
        V=ext.array()
        again=V[:,:idx].copy()
        again[[0,5]]*=np.array([b,b**3])[:,None]
        np.testing.assert_allclose(again,old,rtol=9e-14,atol=0.)
        N,K,E=mom(S)
        e=abs((N.sum()+V[5,idx:].sum())*z**3/target-1.)
        assert e<2e-13
        flags=C["audit"]["solver"]
        assert flags["used_history"]
        H=np.sqrt((TabulatedRadiation(eos,S["T_visible_GeV"]).state(1.)[1]+E.sum()+ext.energy(1.))/3)/MPL
        # The first q64 anchor before these seven blocks is used only to
        # state the momentum cut. This remains a local, bare inverse fold.
        rate=high(S,1./b,12)[1]/(H*E.sum())
        R=dict(T_GeV=S["T_visible_GeV"],snapshot_sha256=prev,
               independent_count_defect=float(e),kernel_energy_defect=C["operator_summary"]["energy_ledger"],
               iterations=flags["newton_iterations"],frozen_bare_high_inverse_energy_over_rho_H=rate,
               peak=metrics([S["entropy_variable_scalar"],S["entropy_variable_fermion"]]))
        if archive:
            Q=archive[i]["snapshot"]
            assert archive[i]["snapshot_sha256"]==digest(Q)
            assert abs(Q["T_visible_GeV"]/S["T_visible_GeV"]-1)<3e-14
            q=np.asarray(Q["lab_momentum_edges_GeV"])
            np.testing.assert_allclose(q,S["lab_momentum_edges_GeV"][:len(q)],rtol=8e-14,atol=0.)
            r= len(q)-1
            QN,QK,QE=mom(Q)
            R.update(shared_count_L1_relative_q32=(np.sum(abs(N[:,:r]-QN),axis=1)/QN.sum(axis=1)).tolist(),
                     q64_added_region_fraction=(N[:,r:].sum(axis=1)/N.sum(axis=1)).tolist())
        rr.append(R)
    return dict(scope="Accepted selected-operator wide continuation; no future q32 shape used as "
                      "a seed in these seven blocks. Frozen high inverse is a local diagnostic, "
                      "not actual high-heavy number or late abundance.",rows=rr)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path);p.add_argument("seed",type=Path)
    p.add_argument("tapes",nargs="+",type=Path)
    p.add_argument("--narrow",type=Path)
    x=p.parse_args()
    print(json.dumps(run(x.eos,x.seed,x.tapes,x.narrow),indent=2,allow_nan=False))
