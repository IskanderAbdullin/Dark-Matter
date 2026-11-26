"""Sensitivity of empty final fermion points to virtual radial panels."""

import argparse
import json
from pathlib import Path

import numpy as np

from joint_visible_heavy_adiabatic_reference import snapshot_functions
from warm_light_occupation_stage_audit_reference import digest,read
from warm_pair_common_fermion_chart_reference import rate


def inspect(path):
    row=read(path); row=row[-1] if isinstance(row,list) else row
    S=row["snapshot"]
    assert row["snapshot_sha256"]==digest(S)
    p=np.asarray(S["lab_momentum_edges_GeV"])
    B,_,F,_=snapshot_functions(S)
    results=[]
    states=[]
    for n in (2,4,8):
        a,g,v,defect,cap=rate(p,F,B,(3,7,6,2),ext=n)
        states.append(g[0])
        results.append(dict(virtual_panels=n,maximum_final_fermion_p=cap,
                            maximum_material_p=float(p[-1]),
                            reverse_rate_number_energy=v[0].tolist(),
                            fraction_inverse_light_products_exterior=float(g[0,1,0,-1]/(2*v[0,0])),
                            ledger=defect))
    grid_L1=(np.abs(states[2]-states[1]).sum(axis=-1)/states[2].sum(axis=-1)).tolist()
    return dict(snapshot_sha256=digest(S), cases=results, four_to_eight_virtual_inverse_grid_L1=grid_L1,
                scope="Only the empty final fermion virtual radial mesh is varied; represented input phase space is fixed.")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("path",type=Path)
    print(json.dumps(inspect(p.parse_args().path),indent=2,allow_nan=False))
