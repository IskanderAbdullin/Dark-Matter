"""Identify native coarse panels that cross target scalar population edges."""
import argparse
import json
from pathlib import Path

import numpy as np

from early_pole_endpoint_cold_normalized_profile_envelope import fit_state,geometry
from early_pole_endpoint_continuum_cold_band_pilot import refine_panels


def build(native,target):
    missing=target[(native[0]<target)&(target<native[-1])&~np.isin(target,native)]
    src=np.searchsorted(native,missing)-1
    if len(src)==0 or len(src)!=len(np.unique(src)):
        raise ArithmeticError('Expected exactly one new target knot per crossing native panel')
    plan=np.c_[native[src],missing,native[src+1]]
    refined,new=refine_panels(native,plan)
    actual=target[(native[0]<target)&(target<native[-1])]
    if not np.all(np.isin(actual,refined)) or len(new)!=2*len(src):
        raise ArithmeticError('Failed to align target occupation with the refined panel set')
    return dict(targeted_parent_panel_splits=plan.tolist(),native_panel_count=len(native)-1,
                refined_panel_count=len(refined)-1,new_coefficients_per_mesh=len(new),
                scope='New coefficients must be directly integrated: plan aligns occupation edges only.')


def run(s8,s16,states,archive):
    x,y=[fit_state(s,archive,n,f) for s,n,f in zip(states,(8,16),(s8,s16),strict=True)]
    plan=[build(geometry(f)[0],y[0].ds.ed) for f in s8]
    if plan[0]!=plan[1]:
        raise ArithmeticError('The two original S8 meshes disagree on target crossings')
    return dict(levels=[8,16],states_sha256_uncompressed=[x[-1],y[-1]],**plan[1])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('s8','s16','states'):
        p.add_argument('--'+k,type=Path,nargs=2,required=True)
    p.add_argument('--archive',type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.s8,a.s16,a.states,a.archive),indent=2,allow_nan=False))
