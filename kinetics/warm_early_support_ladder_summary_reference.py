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

def extend(left,right,finest):
    D=run(left,right)
    F=read(finest)
    assert F["complete"] and F["coarse_and_refined_steps"]==[16,32]
    for row, f in zip(D["matched"],F["matched"][3::4],strict=True):
        assert abs(row["b"]/f["b"]-1.)<1e-14 and abs(row["visible_T_GeV"]/f["T_GeV"]-1.)<1e-12
        L=np.asarray(row["log_support_differences_4_8_16"])
        s16,s32=np.asarray(f["log_support_shift_16_and_32"])
        assert max(abs(L[-1]-s16))<2e-14
        inc0,inc1=np.diff(L,axis=0)
        inc2=s32-s16
        row["log_support_differences_4_8_16_32"]=np.vstack([L,s32]).tolist()
        row["same_support_sign_all_four"]=np.all(L[0]*np.vstack([L,s32])>0,axis=0).tolist()
        row["absolute_last_two_increment_ratios"]=[abs(inc1/inc0).tolist(),abs(inc2/inc1).tolist()]
        row["apparent_last_two_orders"]=[np.log2(abs(inc0/inc1)).tolist(),np.log2(abs(inc1/inc2)).tolist()]
        row["remaining_increment_to_finest_shift"]=abs(inc2/s32).tolist()
    D["inputs"].append(sha(finest))
    D["scope"]="Four completed, exactly shared warm time grids. Signs and empirical differences are numerical diagnostics only, not an extrapolation, late abundance or full support result."
    return D


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("left",type=Path)
    p.add_argument("right",type=Path)
    p.add_argument("--finest",type=Path)
    a=p.parse_args()
    result=extend(a.left,a.right,a.finest) if a.finest else run(a.left,a.right)
    print(json.dumps(result,indent=2,allow_nan=False))
