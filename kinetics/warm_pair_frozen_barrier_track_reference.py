"""Apply the all-tagged-momenta frozen scalar-pair sign bound across saved states."""

import argparse
import json
from pathlib import Path

from warm_light_occupation_stage_audit_reference import digest,read
from warm_pair_frozen_continuum_cap_bound_reference import bound


def track(seed,work,count,style):
    initial=read(seed)
    base=initial[-1] if isinstance(initial,list) else initial
    at=base["snapshot_sha256"]
    assert at==digest(base["snapshot"])
    rows=[bound(seed)]
    for k in range(1,count+1):
        p=work/f"step-{k}.json"
        v=read(p)
        assert v["step" if style=="guided" else "index"]==k
        assert v["parent_snapshot_sha256" if style=="guided" else "previous"]==at
        a=bound(p)
        at=a["snapshot_sha256"]
        rows.append(a)
    return dict(scope="Pointwise across tag momenta on the listed frozen reconstructions, not a trajectory-wide dynamical cap.",
                states_checked=len(rows),minimum_represented_entropy_gap=min(x["entropy_gap_inside"] for x in rows),
                worst_exterior_relative_bound=max(x["exterior_upper_over_interior_lower"] for x in rows),
                uniform_least_negative_net_upper_GeV=max(x["net_upper_GeV"] for x in rows),
                frames=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("seed",type=Path)
    p.add_argument("work",type=Path)
    p.add_argument("count",type=int)
    p.add_argument("--style",choices=("own","guided"),default="own")
    a=p.parse_args()
    print(json.dumps(track(a.seed,a.work,a.count,a.style),indent=2,allow_nan=False))
