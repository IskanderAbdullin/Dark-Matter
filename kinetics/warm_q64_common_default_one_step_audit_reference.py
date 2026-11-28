"""Verify which locked choices reproduce the two published adaptive histories."""

import argparse
import json
from pathlib import Path

import numpy as np

from warm_history_lower_matched_audit_reference import measures
from warm_light_occupation_stage_audit_reference import digest,read


def test(default,one,fine,two):
    A,B,C,D=map(read,(default,one,fine,two))
    rows=[A,B,C,D]
    states=[x["snapshot"] for x in rows]
    for x,s in zip(rows,states,strict=True): assert x["snapshot_sha256"]==digest(s)
    assert B["source"]["guide_sha256"]==A["snapshot_sha256"]
    assert D["source"]["guide_sha256"]==C["snapshot_sha256"]
    assert len(B["logged_stiff_substeps"])==1 and len(D["logged_stiff_substeps"])==2
    assert B["source"]["selected_substeps"]==1 and D["source"]["selected_substeps"]==2
    assert B["source"]["parent_sha256"]==D["source"]["parent_sha256"]
    assert B["source"]["eos_sha256"]==D["source"]["eos_sha256"]
    assert abs(B["logged_stiff_substeps"][0]["step_H"]-sum(z["step_H"] for z in D["logged_stiff_substeps"]))<1e-15
    for s in states:
        np.testing.assert_allclose(s["lab_momentum_edges_GeV"],states[0]["lab_momentum_edges_GeV"],rtol=1e-14)
        np.testing.assert_allclose(s["T_visible_GeV"],states[0]["T_visible_GeV"],rtol=1e-14)
    q=[measures(s) for s in states]
    errors=[]
    for k in (0,2):
        e=abs(q[k]-q[k+1]).sum(axis=-1)/q[k].sum(axis=-1)
        assert np.max(e)<3e-14
        errors.append(e.tolist())
    return dict(scope="Observe that locked one/default and locked two/fine reproduce previous "
                      "adaptive end states to closure precision; not a reconstruction of hidden solver logs.",
                matches=["default one vs adaptive","finer two vs adaptive"],
                step_H=[[s["step_H"] for s in x["logged_stiff_substeps"]] for x in (B,D)],
                normalized_shell_moment_L1=errors,
                states_sha256=[x["snapshot_sha256"] for x in rows])


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("default","locked_one","fine","locked_two"):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    print(json.dumps(test(a.default,a.locked_one,a.fine,a.locked_two),indent=2,allow_nan=False))
