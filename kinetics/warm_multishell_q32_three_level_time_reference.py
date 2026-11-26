"""Compare signed shell errors on 12/24/48 corrected q32 trajectories."""

import argparse
import json
from pathlib import Path

import numpy as np

from warm_multishell_q32_time_audit_reference import series


def run(eos,seed,folders,count):
    q=[]
    m=[]
    for j,D in enumerate(folders):
        record,moment=series(eos,seed,D,12*2**j,count*2**j)
        m.append(moment[::2**j])
        q.append(record)
    for j in (1,2):
        np.testing.assert_allclose(q[0]["T_GeV"],q[j]["T_GeV"][::2**j],rtol=1.5e-13)
    x,y,z=m
    a=x-y
    b=y-z
    N=x[0].sum(axis=-1)
    coarse=abs(a).sum(axis=-1)/N
    fine=abs(b).sum(axis=-1)/N
    p=np.log2(np.divide(coarse[1:],fine[1:]))
    ip=(a[1:]*b[1:]).sum(axis=-1)
    cosine=ip/np.sqrt((a[1:]**2).sum(axis=-1)*(b[1:]**2).sum(axis=-1))
    assert np.all(np.isfinite(p)) and np.max(abs(cosine))<=1.+2e-14
    return dict(scope="Three-grid observed differences; apparent ratios need not be an asymptotic time order.",
                metric_order=["number","full_energy","pressure"],
                shared_coarse_indices=list(range(count+1)),trajectories=q,
                first_level_shell_L1=coarse.tolist(),second_level_shell_L1=fine.tolist(),
                apparent_orders_no_initial=p.tolist(),signed_cosine_no_initial=cosine.tolist())


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for x in ("eos","seed","coarse","medium","fine"):
        p.add_argument(x,type=Path)
    p.add_argument("--limit",type=int,default=12)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.seed,(a.coarse,a.medium,a.fine),a.limit),indent=2,allow_nan=False))
