"""Validate a narrow-solution initial guess without changing the wider target."""

import argparse
import json
from pathlib import Path

import numpy as np

from joint_resonant_centered_stiff_reference import advance
from joint_resonant_frw_split_reference import pack,redshift
from scalar_pair_expanding_continuation_reference import e_total,moments,rate,reconstruct
from warm_coupled_newton_condition_probe_reference import prepare,probe
from warm_light_occupation_stage_audit_reference import digest,read


def run(eos,D,previous,latest):
    a=read(D/"warm-light-occupation-next-q32-tape-2026-11-26.json.gz")
    x,y=read(previous),read(latest)
    S,R=x["snapshot"],y["snapshot"]
    Z=a[2]["snapshot"]
    for o in (x,y,a[2]):
        assert o["snapshot_sha256"]==digest(o["snapshot"])
    assert x["blocks_finished"]==1 and y["blocks_finished"]==2 and x["plan"]==y["plan"]
    b=Z["scale_factor"]/S["scale_factor"]
    assert abs(b/(R["scale_factor"]/S["scale_factor"])-1)<4e-14
    model,u,ref,R0,h,T,H=prepare(eos,S,b)
    _,_,td=model.reference(u,T)
    start=[a.copy() for a in model.arrays(u)]
    size=len(Z["lab_momentum_edges_GeV"])-1
    for i,j in enumerate(("scalar","fermion")):
        start[i][:,:size]=Z["entropy_variable_"+j]
    guess=pack(start)
    history=[]
    def observer(i,r,j):history.append([i,r,j])
    # Imported q32 shape is used *only as an initial guess*. Both input
    # moments and the operator defining the target are exactly q64.
    V,log=advance(model,u,h,td,reference_guess=False,vacuum_scaled=True,
                  initial_theta=guess,balanced_merit=True,observer=observer)
    grid=np.asarray(S["lab_momentum_edges_GeV"])
    p=grid/np.sqrt(b)
    th=model.arrays(V)
    base=moments(p,th)
    pre,_,_=rate(p,th,(3,7,6,2))
    trial=reconstruct(p,base+h*pre[:,:,1:-1]/(4*H))
    signed,_,_=rate(p,trial,(3,7,6,2))
    th=reconstruct(p,base+h*signed[:,:,1:-1]/(2*H))
    th,*_=redshift(grid,np.sqrt(b),b,th)
    U=moments(grid/b,th)
    ref_m=moments(np.asarray(R["lab_momentum_edges_GeV"]),
                  [R["entropy_variable_scalar"],R["entropy_variable_fermion"]])
    size_all=np.array([np.sum(abs(ref_m[:,0]),axis=1),np.full(2,e_total(ref_m))]).T
    dif=np.sum(abs(U-ref_m),axis=2)/size_all
    assert np.max(dif)<4e-10
    assert max(abs(np.sum(U[:,0],axis=1)/np.sum(ref_m[:,0],axis=1)-1.))<1e-11
    return dict(scope="Actual first q64 checkpoint and an independently archived q32 end shape as "
                      "a borrowed initial guess for a single replayed coupled block. The wider "
                      "physical equations, target and outside modes are not replaced. Only the "
                      "retained endpoint is compared here; no production solver changes.",
                wider_initial_sha256=x["snapshot_sha256"],wider_reference_endpoint_sha256=y["snapshot_sha256"],
                narrower_initial_guess_sha256=a[2]["snapshot_sha256"],
                baseline_initial_newton_systems=[probe(model,u,ref,R0,h,i) for i in (False,True)],
                transplanted_iterations_and_residuals=history,solve=log,
                retained_L1_number_and_kinetic_differences_vs_archived=dif.tolist(),
                full_retained_energy_relative_difference=abs(e_total(U)/e_total(ref_m)-1.))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for label in ("eos","previous","latest"):p.add_argument(label,type=Path)
    p.add_argument("--tapes",type=Path,default=Path("results"))
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.tapes,a.previous,a.latest),indent=2,allow_nan=False))
