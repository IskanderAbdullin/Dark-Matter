"""Audit physical endpoint moments on the two-, four- and eight-block q64 ladder."""

import argparse
import json
from pathlib import Path

import numpy as np

from warm_light_occupation_stage_audit_reference import digest,read
from warm_q64_guided_refinement_reference import end_moment


def run(files):
    a=read(files[0])[-1]
    b,c=map(read,files[1:])
    for k in (a,b,c):
        assert k["snapshot_sha256"]==digest(k["snapshot"])
    arr=[]; s=[]
    grid=np.asarray(a["snapshot"]["lab_momentum_edges_GeV"])
    for n,z in zip((2,4,8),(a,b,c),strict=True):
        t=z["snapshot"]
        assert z["step"]==n
        assert abs(t["scale_factor"]/a["snapshot"]["scale_factor"]-1)<5e-14
        np.testing.assert_allclose(t["lab_momentum_edges_GeV"],grid,rtol=8e-14,atol=0.)
        m=end_moment(t)
        arr.append(m)
        s.append(dict(steps=n,temperature_GeV=t["T_visible_GeV"],snapshot_sha256=z["snapshot_sha256"],
                      endpoint_retained_number=m[0].sum(axis=1).tolist(),
                      final_trial=z["audit"]["trial"]["newton_iterations"]))
    norm=np.vstack((arr[0][0].sum(axis=1),np.repeat(arr[0][1].sum(),2)))
    raw=[(abs(x-y).sum(axis=2)/norm) for x,y in zip(arr[:-1],arr[1:],strict=True)]
    assert np.all((raw[0]/raw[1]>1.9)&(raw[0]/raw[1]<2.1))
    return dict(scope="Same physical q64 endpoint, same original operators. Guiding interpolation affects "
                      "only seeds. Observed finite-step error ratio, not an exact continuum result.",
                field_order=["number","full_energy"],controls=s,L1_changes=[x.tolist() for x in raw],
                consecutive_difference_ratio=(raw[0]/raw[1]).tolist())


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("files",type=Path,nargs=3)
    x=p.parse_args()
    print(json.dumps(run(x.files),indent=2,allow_nan=False))
