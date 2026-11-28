"""Compare completed native CLASS cold-limit controls, not a PSD transfer."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator


def load(out,case,redshift_index,k):
    paths=sorted(p for p in out.glob(f"{case}_*_z{redshift_index}_pk.dat")
                 if p.name[len(case)+1:].split("_",1)[0].isdigit())
    a=np.loadtxt(paths[-1])
    assert a[0,0]<=min(k) and a[-1,0]>=max(k)
    return np.exp(PchipInterpolator(np.log(a[:,0]),np.log(a[:,1]))(np.log(k)))


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("class_output",type=Path);a=p.parse_args()
    out=a.class_output
    test=np.array([.01,.05,.1,.25,.5,1.,2.,2.8])
    k=np.unique(np.r_[test,np.geomspace(.01,2.9,400)])
    j=np.searchsorted(k,test)
    rows=[]
    for index,z in enumerate((0,5),1):
        b=load(out,"cdm_true_sync",index,k)
        sync=load(out,"fermion_diagnostic_sync",index,k)/b
        fine=load(out,"fd_true_sync_fine",index,k)/b
        n1=load(out,"fermion_fd_rms",index,k)/b
        n2=load(out,"fermion_fd_rms_fine",index,k)/b
        assert np.max(abs(n1-n2))>.05
        assert np.max(abs(fine-1.))<.02
        rows.append(dict(redshift=z,k_h_per_Mpc=test.tolist(),
            synchronous_standard_over_reference=sync[j].tolist(),
            synchronous_fine_over_reference=fine[j].tolist(),
            Newtonian_standard_over_reference=n1[j].tolist(),Newtonian_fine_over_reference=n2[j].tolist(),
            max_Newtonian_precision_disagreement=float(max(abs(n1-n2))),
            max_synchronous_precision_disagreement=float(max(abs(sync-fine))),
            max_fine_synchronous_deviation_from_CDM=float(max(abs(fine-1.)))))
    print(json.dumps(dict(rows=rows,scope="RMS-momentum-matched analytic FD control in CLASS v3.3.4; numerical failure gate only. Synchronous control includes tiny CDM anchor; not a transfer result for the actual archived spectrum."),indent=2,allow_nan=False))
