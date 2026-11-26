"""Separate time-evolution and integration-sampling shifts in inverse exposure."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid, simpson

from warm_light_occupation_stage_audit_reference import read


def run(paths):
    S=[read(p) for p in paths]
    arrays=[]
    for n,s in zip((4,8,16),S,strict=True):
        assert s["samples_per_decade"]==n
        assert s["source_sha256"]==S[0]["source_sha256"]
        assert s["hard"]==S[0]["hard"] and s["seed"]==S[0]["seed"]
        X,Y=np.asarray(s["H_times"]),np.asarray(s["high_rate_rhoH"])
        assert len(X)==5*n+2 and X[0]==0. and np.all(np.diff(X)>0.) and np.all(Y>=0.)
        assert X[-1]==1e-5 and Y.shape==(len(X),2)
        dose=cumulative_trapezoid(Y,X,axis=0,initial=0.)
        np.testing.assert_allclose(dose,s["integrated_energy_fraction"],rtol=4e-15)
        j=int(np.argmin(abs(X-1e-7)))
        assert abs(X[j]/1e-7-1.)<1e-14
        np.testing.assert_allclose(dose[j]/dose[-1],s["exposure_by_1e_minus_7_over_final"],rtol=4e-15)
        assert np.max(np.abs([r["moment_error"] for r in s["residual_ledgers"]]))<4e-12
        # A quadrature-only estimate from the same solved path. The
        # unsampled first interval is still approximated by its trapezoid.
        log=np.log(X[1:])
        smooth=dose[1]+simpson(Y[1:]*X[1:,None],x=log,axis=0)
        arrays.append((X,Y,dose[-1],dose[j]/dose[-1],smooth))
    comp=[]
    for i in (0,1):
        x,y,I,*_=arrays[i]
        X,Y,J,*_=arrays[i+1]
        ix=np.r_[0,np.arange(1,len(X),2)]
        np.testing.assert_allclose(x,X[ix],rtol=4e-15,atol=0.)
        same_samples=cumulative_trapezoid(Y[ix],x,axis=0)[-1]
        comp.append(dict(coarse_and_fine_steps_per_decade=[4*2**i,8*2**i],
                         total_shift_over_coarse=((J-I)/I).tolist(),
                         sampling_shift_on_finer_path_over_coarse=((J-same_samples)/I).tolist(),
                         path_shift_on_coarse_sampling_over_coarse=((same_samples-I)/I).tolist(),
                         terminal_rate_relative_shift=(Y[-1]/y[-1]-1.).tolist()))
    rows=[dict(samples_per_decade=n,cumulative_energy_fraction=A[2].tolist(),
               fraction_through_Hdt_1e_minus_7=A[3].tolist(),terminal_normalized_rate=A[1][-1].tolist(),
               log_simpson_same_path_fraction=A[4].tolist()) for n,A in zip((4,8,16),arrays,strict=True)]
    return dict(scope="Same isolated fixed clock and hard-packet seed. Numerical quadrature "
                      "decomposition, not a comoving dose or rigorous stronger bound.",
                hard_GeV=S[0]["hard"],source_snapshot_sha256=S[0]["source_sha256"],
                rate_order=["bare_inverse_energy","tagged_inverse_energy"],ladders=rows,
                decompositions=comp)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("per_decade4","per_decade8","per_decade16"):
        p.add_argument(name,type=Path)
    a=p.parse_args()
    print(json.dumps(run((a.per_decade4,a.per_decade8,a.per_decade16)),indent=2,allow_nan=False))
