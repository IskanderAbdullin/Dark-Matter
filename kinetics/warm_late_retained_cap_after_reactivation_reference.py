"""Audit whether the original retained-shell ceiling persists after activation."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import expit, roots_legendre

from joint_prompt_two_moment_evolution_reference import kin
from three_species_decay_onsager_reference import m,mp
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom


def f(a,g):
    return expit(-a) if g==4 else np.exp(-a)/(-np.expm1(-a))


def run(seed,paths):
    src=read(seed)
    root=src["snapshot"]
    assert src["snapshot_sha256"]==digest(root)
    p=np.asarray(root["lab_momentum_edges_GeV"])
    n=len(p)-1
    keys=["entropy_variable_"+x for x in ("scalar","fermion")]
    original=[np.asarray(root[k]) for k in keys]
    N=get_mom(p,original)[0].sum(axis=1)
    z,w=roots_legendre(18)
    out=[]
    for path in paths:
        a=read(path)
        S=a["snapshot"]
        assert a["hash"]==digest(S)
        pp=np.asarray(S["lab_momentum_edges_GeV"])
        np.testing.assert_array_equal(p,pp[:n+1])
        counts=[]
        for j,(key,mass,g) in enumerate(zip(keys,(mp,m),(1,4),strict=True)):
            O=original[j]
            X=np.asarray(S[key])[:,:n]
            cap=O[0]+np.minimum(0.,O[1])
            a0,b=X
            gain=cap-(a0+np.minimum(0.,b))
            eligible=np.flatnonzero(gain>2e-13)
            k=kin(p,mass)
            low=k[:-1].copy()
            high=k[1:].copy()
            q=np.divide(cap-a0,b,out=np.zeros(n),where=b!=0.)
            split=low+np.diff(k)*np.clip(q,0.,1.)
            lower=np.where(b<0.,split,low)[eligible]
            upper=np.where(b>0.,split,high)[eligible]
            Xk=lower[:,None]+(upper-lower)[:,None]*(1+z)/2
            weights=(upper-lower)[:,None]*w/2
            u=(Xk-k[eligible,None])/np.diff(k)[eligible,None]
            flux=f(a0[eligible,None]+b[eligible,None]*u,g)-f(cap[eligible,None],g)
            assert np.min(flux,initial=0.)>=-1e-17
            phi=np.sqrt(Xk*(Xk+2*mass))
            excess=g/(2*np.pi**2)*np.sum(weights*(Xk+mass)*phi*np.maximum(flux,0.))
            rel=np.divide(f(a0+np.minimum(0.,b),g),f(cap,g),out=np.full(n,np.nan),where=f(cap,g)>0.)
            counts.append(dict(fraction_of_initial_retained_count_above_ceiling=float(excess/N[j]),
                               cells_exceeding_original_envelope=len(eligible),
                               biggest_occupied_cell_maximum_over_original=float(np.nanmax(rel)),
                               largest_reduction_in_minimum_entropy=float(np.max(gain))))
        out.append(dict(source_sha256=a["hash"],H_times_elapsed=a["time"],species=counts))
    return dict(scope="Counts only the occupation above the original pointwise ceiling in retained cells. "
                      "Positive violation means the instantaneous anchor upper cannot be presumed persistent.",
                original_snapshot_sha256=src["snapshot_sha256"],selected_states=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("seed",type=Path)
    p.add_argument("paths",nargs="+",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.seed,a.paths),indent=2,allow_nan=False))
