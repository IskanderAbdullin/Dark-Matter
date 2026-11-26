"""Read-only audit of explicit states in one new serial warm block."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.special import expit

import coupled_warm_pair_continuation_reference as kernel
import warm_coupled_serial_restart_reference as serial


def metrics(state):
    a,b=np.asarray(state[0])
    c,d=np.asarray(state[1])
    z=min(np.min(a),np.min(a+b))
    zf=min(np.min(c),np.min(c+d))
    assert z>0
    return dict(retained_scalar_peak=float(np.exp(-z)/(-np.expm1(-z))),
                retained_light_peak=float(expit(-zf)),
                entropy_array_sha256=hashlib.sha256(serial.bytes_of(
                    [np.asarray(t).tolist() for t in state])).hexdigest())


def observe(eos,source,directory,stop,count,outfile):
    tape=[]
    rd,rc,cl=kernel.redshift,kernel.reconstruct,kernel.collide
    seq={"redshift":0,"kick":0}

    def redshift(q,b1,b2,state,*args,**kwargs):
        seq["redshift"]+=1
        z=seq["redshift"]
        tape.append(dict(operation=f"redshift_{z}_input",energy_grid_sha256=hashlib.sha256(
                         serial.bytes_of((np.asarray(q)/b1).tolist())).hexdigest(),**metrics(state)))
        value=rd(q,b1,b2,state,*args,**kwargs)
        tape.append(dict(operation=f"redshift_{z}_output",energy_grid_sha256=hashlib.sha256(
                         serial.bytes_of((np.asarray(q)/b2).tolist())).hexdigest(),**metrics(value[0])))
        return value

    def reconstruct(q,value,*args,**kwargs):
        seq["kick"]+=1
        state=rc(q,value,*args,**kwargs)
        tape.append(dict(operation=f"pair_reconstruction_{seq['kick']}",
                         trial_midpoint=bool(seq["kick"]%2),**metrics(state)))
        return state

    def collide(model,u,h,T,*args,**kwargs):
        value=cl(model,u,h,T,*args,**kwargs)
        th=model.arrays(value[0])
        tape.append(dict(operation="implicit_collision_output",
                         microsteps_accepted=int(value[2]),microsteps_rejected=int(value[3]),
                         **metrics(th)))
        return value

    kernel.redshift,kernel.reconstruct,kernel.collide=redshift,reconstruct,collide
    try:
        serial.main(eos,source,directory,stop,count,True,1)
    finally:
        kernel.redshift,kernel.reconstruct,kernel.collide=rd,rc,cl
    assert seq=={"redshift":2,"kick":4} and len(tape)==9
    ck=json.loads((directory/"checkpoint.json").read_bytes())
    assert tape[-1]["entropy_array_sha256"]==metrics([ck["snapshot"]["entropy_variable_scalar"],
                 ck["snapshot"]["entropy_variable_fermion"]])["entropy_array_sha256"]
    serial.atomic(outfile,dict(blocks_finished=ck["blocks_finished"],plan_sha256=ck["plan_sha256"],
                               checkpoint_sha256=serial.sha(directory/"checkpoint.json"),samples=tape))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for a in ("eos","source","directory","observer_result"):
        p.add_argument(a,type=Path)
    p.add_argument("--stop",type=float,default=1.01)
    p.add_argument("--blocks",type=int,default=2)
    a=p.parse_args()
    observe(a.eos,a.source,a.directory,a.stop,a.blocks,a.observer_result)
