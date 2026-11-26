"""Check that a strict strip is fully masked for one clock window."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from disjoint_hard_thermal_measured_source_reference import WINDOWS


def run(older,later):
    A,B=[json.loads(p.read_text()) for p in (older,later)]
    assert A["mass_over_T_interval"]==[4.,8.] and B["mass_over_T_interval"]==[4.,9.]
    rows=[]
    for a,b in zip(A["controls"],B["controls"]):
        T0=a["Born_T_min_GeV"]
        T1=a["measured_plus_disjoint_Born_root_GeV"]
        assert T0==b["Born_T_min_GeV"] and T1==b["measured_plus_disjoint_Born_root_GeV"]
        # All new-strip masses are inside the first window exactly iff
        # both opposite rectangle corners lie inside that window.
        mn,mx=8*T0,9*T1
        fully=bool(WINDOWS[0,0]<=mn<=mx<=WINDOWS[0,1])
        checks=[]
        for x,y in zip(a["results"],b["results"]):
            mode=x["interpolation"]
            assert mode==y["interpolation"]
            changes=[]
            for p,q in zip(x["alpha_cases"],y["alpha_cases"]):
                assert p["alpha_s"]==q["alpha_s"]
                changes.append(dict(alpha_s=p["alpha_s"],root_difference_GeV=q["T_R_GeV"]-p["T_R_GeV"],
                    masked_increment_through_8=p["fraction_of_strict_increment_cut_for_overlap"],
                    masked_increment_through_9=q["fraction_of_strict_increment_cut_for_overlap"]))
            if fully and mode=="linear":
                assert np.max([abs(c["root_difference_GeV"]) for c in changes])<1e-11
            checks.append(dict(interpolation=mode,changes=changes))
        rows.append(dict(lower_T_GeV=T0,Born_root_upper_T_GeV=T1,new_mass_range_GeV=[mn,mx],
                    completely_inside_lower_measured_window=fully,conditional_root_checks=checks))
    return dict(first_measured_window_GeV=WINDOWS[0].tolist(),mass_over_T_new_strip=[8.,9.],cases=rows,
        result_sha256=[hashlib.sha256(k.read_bytes()).hexdigest() for k in (older,later)],
        scope="Geometric masking of the new strict strip in a bounded optimistic clock integral, not convergence of full thermal QCD. Piecewise linear interpolation is unchanged below x=8; PCHIP can re-slope adjacent intervals even when the extension itself is masked.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("previous",type=Path)
    p.add_argument("extended",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.previous,a.extended),indent=2,allow_nan=False))
