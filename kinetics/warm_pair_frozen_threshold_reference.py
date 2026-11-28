"""Support-uniform tagged occupation threshold from certified scalar-pair rates."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import expit

from warm_history_lower_matched_audit_reference import measures
from warm_light_occupation_stage_audit_reference import digest,read
from warm_pair_frozen_continuum_cap_bound_reference import bound


def threshold(data):
    g=data.get("entropy_gap_inside",data.get("represented_gap"))
    # L_min(1/2) = (1-exp(-g)) times a lower bound on the
    # represented-partner forward rate at occupancy one half.
    ratio=data.get("exterior_upper_over_interior_lower",data.get("relative_upper"))
    assert g>0. and ratio>=0.
    formal=np.exp(-g)+ratio*(-np.expm1(-g))
    return float(formal/(1.+formal))


def stages(path):
    t=read(path)
    out={}
    for key,row in t.items():
        out[key]={}
        for name,v in (("frozen",row["frozen"]["frames"]),("interpolant",row["interpolant"]["brackets"])):
            gt=[threshold(x) for x in v]
            out[key][name]=dict(guaranteed_maximum_occupation_threshold=max(gt),thresholds=gt)
    return out


def archive(path):
    row=read(path)
    row=row[-1] if isinstance(row,list) else row
    S=row["snapshot"]
    assert digest(S)==row["snapshot_sha256"]
    b=bound(path)
    f=threshold(b)
    x=np.asarray(S["entropy_variable_fermion"])
    minimum=expit(-np.maximum(x[0],x[0]+x[1]))
    maximum=expit(-np.minimum(x[0],x[0]+x[1]))
    n=measures(S)[0,1]
    keep=minimum>=f
    any_bin=maximum>=f
    return dict(snapshot_sha256=digest(S),T_GeV=S["T_visible_GeV"],uniform_threshold=f,
                total_retained_shells=len(n),entirely_above_threshold_shells=int(np.sum(keep)),
                first_shell_minimum=float(minimum[0]),first_shell_maximum=float(maximum[0]),
                guaranteed_nonpositive_pair_retained_number_fraction_lower=float(n[keep].sum()/n.sum()),
                possible_above_threshold_retained_number_fraction_upper=float(n[any_bin].sum()/n.sum()))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("tracks",type=Path)
    p.add_argument("archives",nargs="*",type=Path)
    a=p.parse_args()
    print(json.dumps(dict(scope="A uniform sufficient sign threshold for the scalar-pair term with the partner distributions frozen.",
                          track_thresholds=stages(a.tracks),direct_archive_checks=[archive(q) for q in a.archives]),
                     indent=2,allow_nan=False))
