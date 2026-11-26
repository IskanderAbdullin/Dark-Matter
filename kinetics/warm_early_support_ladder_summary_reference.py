"""Combine equal-clock 4–8 and 8–16 support differences at shared points."""
import argparse
import json
from pathlib import Path

import numpy as np

from warm_early_support_time_interaction_reference import read,sha


def run(left,right):
    A,B=map(read,(left,right))
    assert A["complete"] and B["complete"] and B["coarse_and_refined_steps"]==[8,16]
    rows=[]
    for a,b in zip(A["matched"],B["matched"][1::2],strict=True):
        assert a["b"]==b["b"] and abs(a["T_GeV"]/b["T_GeV"]-1.)<1e-12
        D4,D8=np.asarray(a["log_support_shift_4_and_8"])
        d8,D16=np.asarray(b["log_support_shift_8_and_16"])
        e=float(max(abs(D8-d8)))
        assert e<2e-14
        x=D8-D4
        y=D16-d8
        rows.append(dict(b=a["b"],visible_T_GeV=a["T_GeV"],log_support_differences_4_8_16=[z.tolist() for z in (D4,D8,D16)],
                        same_support_sign_all=((D4*D8>0.)&(D8*D16>0.)).tolist(),
                        absolute_delta_8_to_16_over_delta_4_to_8=abs(y/x).tolist(),
                        apparent_refinement_order=np.log2(abs(x/y)).tolist(),shared_eight_block_defect=e))
    return dict(matched=rows,inputs=[sha(left),sha(right)],
           scope="Empirical relative time refinement of the support contrast at four shared warm endpoints. No extrapolated continuum, late relic, or inherited exterior momentum distribution.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("left",type=Path)
    p.add_argument("right",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.left,a.right),indent=2,allow_nan=False))
