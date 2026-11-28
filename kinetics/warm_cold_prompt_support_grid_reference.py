"""Compare actual retained prompt geometry to the independent fixed-light kernel."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import roots_legendre

from joint_prompt_two_moment_evolution_reference import Collision
from joint_visible_heavy_adiabatic_reference import parent_branch,snapshot_functions
from scalar_pair_expanding_continuation_reference import moments
from three_species_decay_onsager_reference import m,mp,width
from warm_light_occupation_stage_audit_reference import digest,read
from warm_near_rest_light_inverse_decay_reference import E0,rate,window


def independent_first_average(S,H,n):
    B,be,F,fe=snapshot_functions(S)
    e=S["lab_momentum_edges_GeV"][1]
    z,w=roots_legendre(n)
    p=e*(z+1)/2
    def heavy(E):
        return parent_branch(E.ravel(),B,F,be,fe,angle=21)[:,0].reshape(E.shape)
    C,*_=rate(p,F,B,heavy,be,order=14)
    net=C[0]+C[1]-C[2]
    return float((w*p*p)@net/np.sum(w*p*p)*width/H)


def geometry(S,H,n):
    p=np.asarray(S["lab_momentum_edges_GeV"])
    shapes=[np.asarray(S["entropy_variable_"+a]) for a in ("scalar","fermion")]
    A=Collision(p,n,n)
    touched=np.bincount(A.J,minlength=len(p)-1)
    N=moments(p,shapes)[1,0]
    f=A.rate(shapes)[1,0]/(A.c[1].vol*H)
    cutoff=np.hypot(p[-1],mp)
    if cutoff<E0:
        # Scalar minimum falls through the cut well inside this grid.
        p0=brentq(lambda x:window(x)[-2]-cutoff,0.,p[-1],xtol=5e-16)
    else:
        p0=0.
    unreachable=np.flatnonzero(p[1:]<=p0)
    assert np.all(touched[unreachable]==0)
    # Independent kinematics for classification of fully disconnected bins;
    # any partially touched cell beyond these may be quadrature-dependent.
    assert np.array_equal(np.flatnonzero(touched==0),unreachable)
    return dict(order=n,first_cell_events=int(touched[0]),
                count_of_fully_inactive_light_shells=len(unreachable),
                first_kinematically_allowed_light_p_GeV=p0,
                summed_fully_inactive_light_number_fraction=float(N[unreachable].sum()/N.sum()),
                projected_first_cell_occupation_rate_per_Hubble=float(f[0]))


def run(D):
    ref=read(D/"warm-near-rest-light-inverse-decay-2026-11-26.json.gz")
    out=[]
    for q in (16,32):
        S=read(D/f"warm-light-occupation-next-q{q}-tape-2026-11-26.json.gz")
        S+=[x["checkpoint"] for x in read(D/f"warm-observed-suboperator-q{q}-tape-2026-11-26.json.gz")["blocks"]]
        z=next(y for y in ref["cases"] if y["q"]==q)["stages"]
        rows=[]
        for c,a in zip(S,z,strict=True):
            s=c["snapshot"];h=c["snapshot_sha256"]
            assert h==digest(s)==a["snapshot_sha256"]
            H=a["fine"]["H_over_rest_decay_width"]*width
            g=[geometry(s,H,n) for n in (3,5,8,12)]
            mean=[independent_first_average(s,H,n) for n in (12,24)]
            assert abs(mean[0]/mean[1]-1)<1e-6
            r=g[-1]["projected_first_cell_occupation_rate_per_Hubble"]
            if q==16:
                assert r==0. and mean[-1]>0.
            else:
                assert r<0. and mean[-1]<0.
                assert abs(r/mean[-1]-1)<2e-4
            rows.append(dict(T_GeV=a["T_GeV"],snapshot_sha256=h,orders=g,
                             independently_integrated_open_decay_first_cell_rate_per_Hubble=mean))
        out.append(dict(q=q,stages=rows))
    return dict(scope="Frozen instant prompt geometry. The prompt solver intentionally requires both "
                      "daughters to be represented; the independent check allows spontaneous decay "
                      "into an unrepresented scalar from a retained stationary H surrogate. "
                      "Agreement on q32 in this small region is not a global open-support result.",
                cases=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tapes",type=Path,default=Path("results"))
    print(json.dumps(run(p.parse_args().tapes),indent=2,allow_nan=False))
