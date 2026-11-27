"""Audit the two independent late common-chart outputs against the nine-cell baseline."""

import argparse
import json
from pathlib import Path
import numpy as np

from frw_exterior_restart_reference import import_restart
from three_species_decay_onsager_reference import m,mp
from warm_history_lower_matched_audit_reference import measures
from warm_light_occupation_stage_audit_reference import digest,read


def run(parent,baseline,default,fine):
    A,B,C,D=map(read,(parent,baseline,default,fine))
    rows=[A,B,C,D]
    states=[x["snapshot"] for x in rows]
    assert all(digest(s)==x["snapshot_sha256"] for x,s in zip(rows,states,strict=True))
    S,G=states[:2]
    assert B["parent_snapshot_sha256"]==A["snapshot_sha256"]
    h=C["source"]["whole_chart_sha256"]
    assert C["source"]==dict(D["source"],fine=False)
    assert C["source"]["parent_snapshot_sha256"]==A["snapshot_sha256"]
    assert C["source"]["preconditioner_sha256"]==B["snapshot_sha256"]
    assert C["source"]["eos_sha256"]==B["plan"]["eos_sha256"]
    assert not C["source"]["fine"] and D["source"]["fine"]
    for s in states[2:]:
        np.testing.assert_allclose(s["lab_momentum_edges_GeV"],G["lab_momentum_edges_GeV"],rtol=3e-14)
        np.testing.assert_allclose(s["T_visible_GeV"],G["T_visible_GeV"],rtol=3e-14)
        np.testing.assert_allclose(s["scale_factor"],G["scale_factor"],rtol=3e-14)
    scaled=[r["scale_factor"]/S["scale_factor"] for r in states[1:]]
    src=import_restart(S["positive_exterior_restart"]).array()
    out=[]
    X=[measures(r) for r in states]
    def cold(M,s):
        return M[0,1,0]/(4*s["lab_momentum_edges_GeV"][1]**3/(6*np.pi**2))
    for label,r,s,z,b in zip(("nine-cell","common default","common fine"),rows[1:],states[1:],X[1:],scaled,strict=True):
        N=import_restart(s["positive_exterior_restart"]).array()
        prefix=N[:,:src.shape[1]].copy()
        prefix[[0,5]]*=np.array([b,b**3])[:,None]
        np.testing.assert_allclose(prefix,src,rtol=1.4e-13,atol=0.)
        birth=N[:,src.shape[1]:]
        receipt=np.array([birth[5,birth[1]==d].sum() for d in (mp,m)])
        out.append(dict(label=label,snapshot_sha256=r["snapshot_sha256"],
                        cold_slope_per_ln_a=float((cold(z,s)-cold(X[0],S))/np.log(b)),
                        new_positive_outside_count_over_retained=(receipt/z[0].sum(axis=1)).tolist(),
                        new_positive_receipt_count=birth.shape[1]))
    norm=X[1].sum(axis=-1)
    return dict(scope="Actual matched one-step solutions. Common-chart runs may use different fallback "
                      "microstep paths, so differences need not isolate quadrature.",
                original_chart_sha256=h,branches=out,
                whole_shell_L1_default_vs_fine=(abs(X[2]-X[3]).sum(axis=-1)/norm).tolist(),
                common_fallback_reasons=[r["solver"] for r in (C,D)])


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("parent","baseline","default","fine"):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.parent,a.baseline,a.default,a.fine),indent=2,allow_nan=False))
