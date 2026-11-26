"""Audit the native cold-light elastic number projection against tagged integration."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from joint_below_pole_elastic_linear_modes_reference import setup
from joint_prompt_two_moment_evolution_reference import Collision
from joint_resonant_frw_split_reference import pack
from joint_visible_heavy_adiabatic_reference import snapshot_functions
from late_archived_shape_continuation_reference import Late
from three_species_decay_onsager_reference import width
from warm_cold_light_below_pole_elastic_reference import collision
from warm_light_occupation_stage_audit_reference import digest,read


def independent(s,H,radial,kernel):
    B,be,F,fe=snapshot_functions(s)
    pe=s["lab_momentum_edges_GeV"][1]
    x,w=roots_legendre(radial)
    p=pe*(x+1)/2
    a=np.array([collision(d,B,be,F,fe,.82,kernel)[0] for d in p])
    result=(w*p*p)@a/np.sum(w*p*p)/H
    assert result[0]<result[1]
    return dict(radial=radial,kernel=kernel,gain_loss_per_Hubble=result.tolist(),
                net_per_Hubble=float(result[0]-result[1]))


def native(s,H,orders):
    grids=np.asarray(s["lab_momentum_edges_GeV"])
    st=pack([s["entropy_variable_scalar"],s["entropy_variable_fermion"]])
    seed=Late(Collision(grids),H,[],True)
    ref,_,T=seed.reference(st,s["T_visible_GeV"])
    cloud=[]
    setup(max(T,s["T_visible_GeV"]),seed.k,T,seed.eta,*orders,32.,.82,
          cache=cloud,cache_only=True)
    op=Late(seed.k,H,cloud,True,compact=True)
    R,_,T2=op.reference(st,s["T_visible_GeV"])
    assert np.array_equal(ref,R) and T==T2
    d=op.elastic(st-R)[0].reshape(2,seed.n,2)
    N=d[:,:,0]
    energy=np.array([N[j]@(c.mass+c.low)+d[j,:,1]@c.dk for j,c in enumerate(seed.k.c)])
    den=np.sum(abs(N),axis=1)
    defect=np.r_[abs(N.sum(axis=1))/den,abs(energy.sum())/np.sum(abs(energy))]
    assert max(defect)<5e-11
    return dict(orders=orders,panels=len(cloud),
                first_light_average_per_Hubble=float(N[1,0]/seed.k.c[1].vol[0]),
                conserved_number_and_energy_residuals=defect.tolist())


def run(D):
    R=read(D/"warm-near-rest-light-inverse-decay-2026-11-26.json.gz")
    rows=[]
    for q in (16,32):
        t=read(D/f"warm-light-occupation-next-q{q}-tape-2026-11-26.json.gz")
        a=next(z["stages"] for z in R["cases"] if z["q"]==q)
        for j in (0,3):
            s=t[j]["snapshot"]
            assert digest(s)==t[j]["snapshot_sha256"]==a[j]["snapshot_sha256"]
            H=a[j]["fine"]["H_over_rest_decay_width"]*width
            I=[independent(s,H,p,op) for p,op in ((3,(6,6,6,6)),(6,(9,8,9,8)))]
            N=[native(s,H,op) for op in ((3,2,2,8),(4,3,3,10),(6,4,4,12),(9,6,6,16))]
            assert abs(I[1]["net_per_Hubble"]/I[0]["net_per_Hubble"]-1)<.004
            assert max(v["first_light_average_per_Hubble"] for v in N)<0.
            rows.append(dict(q=q,j=j,T_GeV=s["T_visible_GeV"],snapshot_sha256=t[j]["snapshot_sha256"],
                             independent=I,native=N,
                             native_fine_over_independent_fine=N[-1]["first_light_average_per_Hubble"]/I[-1]["net_per_Hubble"]))
    return dict(scope="Selected below-pole elastic operator on frozen warm states only. "
                      "The continuum tagged integral uses the archived vacuum exterior assumption; "
                      "the native quadrature projects conservative redistributions on retained "
                      "support and imposes its existing energy and thermal panel conventions. "
                      "Agreement in the first cold cell does not certify global grid convergence.",
                cases=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tapes",type=Path,default=Path("results"))
    print(json.dumps(run(p.parse_args().tapes),indent=2,allow_nan=False))
