"""Contrast identical CLASS runs changing only the perturbation derivative."""
import argparse
import json
from pathlib import Path
import re

import numpy as np

from archived_class_pressure_balanced_reference import read,at
from archived_class_manual_momentum_gate_reference import audit


def difference(out,raw,patched):
    check=re.search(r'DILATION_DIAGNOSTIC n=(\d+) raw_residual=(\S+) factor=(\S+)',
                    (out/(patched+".log")).read_text())
    if not check:
        raise ValueError("No diagnostic marker: expected perturbed build")
    n,d,c=check.groups()
    n,d,c=int(n),float(d),float(c)
    assert abs((1+d)*c-1.)<2e-13
    bg0=np.loadtxt(out/(raw+"_00_background.dat"))
    bg1=np.loadtxt(out/(patched+"_00_background.dat"))
    assert np.array_equal(bg0,bg1)
    rows=[]
    for i,z in enumerate((0,5),1):
        A=read(out,raw,i)
        B=read(out,patched,i)
        C=read(out,"cdm_balanced_grid_sync",i)
        assert np.array_equal(A[:,0],B[:,0])
        mask=(C[:,0]>=.012)&(C[:,0]<=2.9)
        k=C[mask,0]
        P0=at(A,k)/C[mask,1]
        P1=at(B,k)/C[mask,1]
        rows.append(dict(z=z,mode_count=len(k),k_bounds_h_per_Mpc=[k[0],k[-1]],
            raw_max_CDM_departure=float(np.max(np.abs(P0-1.))),
            diagnostic_max_CDM_departure=float(np.max(np.abs(P1-1.))),
            max_paired_change=float(np.max(np.abs(P0-P1))),
            grid_h_per_Mpc=k.tolist(),raw_over_CDM=P0.tolist(),
            derivative_rescaled_over_CDM=P1.tolist()))
    return dict(nodes=n,native_number_dilation_residual=d,applied_derivative_factor=c,
                identical_background_array=True,rows=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("out",type=Path)
    p.add_argument("fine_table",type=Path)
    p.add_argument("--nodes",nargs="+",type=int,default=[64,128])
    a=p.parse_args()
    cases=[difference(a.out,f"psd_raw{n}_matched_sync",f"psd_dilation{n}_sync") for n in a.nodes]
    r=audit(a.fine_table,a.nodes,20.)["results"]
    for case,z in zip(cases,r):
        d=z["dilation_integral_by_parts_residual"][0]
        case["independent_cubic_spline_number_dilation_residual"]=d
        assert abs(case["native_number_dilation_residual"]-d)<1e-6
    print(json.dumps(dict(restricted_cases=cases,
        scope="Native-control differencing with a diagnostic-only change to the weighted derivative. Not a physical transfer, replacement CLASS algorithm, upstream correction or validation of the modified solver."),
        indent=2,allow_nan=False))
