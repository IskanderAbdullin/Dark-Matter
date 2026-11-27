"""Audit sign brackets, total discrete transfer and selected mesh/clock changes."""

import argparse
import json
from pathlib import Path

import numpy as np

from warm_light_occupation_stage_audit_reference import read


def metrics(v):
    r=v["rows"]; plan=r[0]["plan"]
    n=plan["locked_microsteps"]
    assert len(r)==n+1 and all(x["plan"]==plan for x in r)
    steps=np.array([x["time"] for x in r])
    h=(steps[-1]-steps[0])/n
    np.testing.assert_allclose(np.diff(steps),h,rtol=3e-14)
    W=np.array([x["tail_energy"] for x in r])
    rate=np.array([x["signed_native_flux"] for x in r])
    P=np.array([x["angular_consumed_and_emitted_daughter_excess"] for x in r])
    d=np.array([x["exact_discrete_tail_balance_residual"] for x in r])
    np.testing.assert_allclose(np.diff(W)-rate[1:]*h,d[1:],rtol=0.,atol=1e-32)
    assert abs(W[-1]-W[0]-h*rate[1:].sum())<2e-10*np.max(W)
    assert all(x["accepted"]["residual"]<8e-12 for x in r[1:])
    bracket=np.flatnonzero(rate[1:]*rate[:-1]<0.)
    G=[list(steps[k:k+2]) for k in bracket]
    assert G==v["sign_change_brackets"]
    return steps,W,rate,dict(substeps=n,sign_brackets=G,max_sampled_excess=float(W.max()),
                             first_maximum_sample_time=float(steps[np.argmax(W)]),
                             signed_clock_integral=float(h*rate[1:].sum()),
                             independent_angular_gross_scale_error=float(np.max(abs(rate-(P[:,1]-P[:,0]))/P.sum(axis=1))),
                             gross_consumed_at_start_and_end=[float(P[0,0]),float(P[-1,0])],
                             gross_emitted_at_start_and_end=[float(P[0,1]),float(P[-1,1])],
                             accumulated_discrete_ledger_residual=float(W[-1]-W[0]-h*rate[1:].sum()))


def run(paths,mesh=None):
    vals=[read(p) for p in paths]
    keys=[x["rows"][0]["plan"]["seed_sha256"] for x in vals]
    assert len(set(keys))==1
    M=[metrics(x) for x in vals]
    controls=[]
    for a,b in zip(M[:-1],M[1:],strict=True):
        assert len(b[0])-1==2*(len(a[0])-1)
        np.testing.assert_allclose(a[0],b[0][::2],rtol=3e-14)
        controls.append(dict(coarse_and_fine_substeps=[len(a[0])-1,len(b[0])-1],
                             max_matched_excess_difference_over_fine_peak=float(np.max(abs(a[1]-b[1][::2]))/b[1].max())))
    out=dict(scope="Exact discrete ledger on the selected locked interval. Sign brackets are sampled "
                   "properties of these time-discrete trajectories, not a continuum enclosure.",
             native_seed_sha256=keys[0],results=[r[3] for r in M],refinements=controls)
    if mesh:
        v=read(mesh)
        t,y,g,D=metrics(v)
        assert len(t)==len(M[0][0])
        np.testing.assert_allclose(t,M[0][0],rtol=3e-14)
        np.testing.assert_allclose(y[0],M[0][1][0],rtol=3e-13)
        out["split_mesh_control"]=dict(**D,
                                        max_matched_excess_difference_over_original_peak=float(np.max(abs(y-M[0][1]))/M[0][1].max()))
    return out


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("ladders",type=Path,nargs="+")
    p.add_argument("--mesh",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.ladders,a.mesh),indent=2,allow_nan=False))
