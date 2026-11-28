"""Prescribed warm grid with atomic positive-cohort checkpoints at outer blocks."""
import argparse
import fcntl
import hashlib
import json
import os
import sys
import tempfile
from pathlib import Path

import numpy as np

from coupled_warm_pair_continuation_reference import run
from joint_nonlinear_prompt_return_moments_reference import P


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def bytes_of(x):
    return (json.dumps(x,sort_keys=True,separators=(",",":"),allow_nan=False)+"\n").encode()


def atomic(path, value):
    # Old and new path are on the same filesystem. Persist before advancing.
    fd,tmp=tempfile.mkstemp(prefix=path.name+".",dir=path.parent)
    try:
        with os.fdopen(fd,"wb") as f:
            f.write(bytes_of(value))
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp,path)
        root=os.open(path.parent,os.O_RDONLY)
        try:
            os.fsync(root)
        finally:
            os.close(root)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def main(eos, start, directory, stop, count, local, max_blocks):
    if count<1 or not np.isfinite(stop) or stop<=1. or (max_blocks is not None and max_blocks<1):
        raise ValueError("Require positive block count and finite expansion greater than 1.")
    base=json.loads(start.read_text())["snapshots"][-1]
    grid=np.unique(np.r_[np.linspace(0.,np.log(stop),count+1),
                         np.clip(np.log(base["lab_momentum_edges_GeV"][-1]/P),0.,np.log(stop))])
    assert np.all(np.diff(grid)>1e-14)
    code=Path(__file__).resolve().parent
    modules=[Path(__file__).resolve(),code/"coupled_warm_pair_continuation_reference.py",
             code/"frw_exterior_restart_reference.py"]
    plan=dict(reference_sha256={p.name:sha(p) for p in modules},schema="warm-fixed-log-grid-positive-cohorts-v1", source_sha256=sha(start), eos_sha256=sha(eos),
              outer_log_scale_grid=grid.tolist(),local_prompt_jacobian=bool(local))
    plan_hash=hashlib.sha256(bytes_of(plan)).hexdigest()
    directory.mkdir(parents=True,exist_ok=True)
    with (directory/"workspace.lock").open("a+b") as lock:
        try:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            raise ValueError("A solver still owns this work directory.") from None
        saved=directory/"checkpoint.json"
        cur=base
        start_index=0
        history=[]
        if saved.exists():
            state=json.loads(saved.read_text())
            if state["plan"]!=plan or state["plan_sha256"]!=plan_hash:
                raise ValueError("Saved checkpoint is for different inputs, grids, or solver flags.")
            history=state["ledger"]
            start_index=state["blocks_finished"]
            assert 0<start_index<len(grid) and start_index==len(history)
            assert state["blocks_total"]==len(grid)-1
            for j,row in enumerate(history,1):
                assert row["block"]==j and abs(row["relative_expansion"]/np.exp(grid[j])-1.)<3e-12
            cur=state["snapshot"]
            assert abs(cur["scale_factor"]/base["scale_factor"]/np.exp(grid[start_index])-1.)<3e-12
            assert state["snapshot_sha256"]==hashlib.sha256(bytes_of(cur)).hexdigest()
        step_end=len(grid)-1 if max_blocks is None else min(len(grid)-1,start_index+max_blocks)
        for j in range(start_index,step_end):
            stop_j=float(np.exp(grid[j+1]-grid[j]))
            nxt={}
            r=run(eos,cur,1,stop_j,checkpoint=nxt,track_exterior=True,local_prompt_jacobian=local,
                  scheduled_single_block=True)
            assert len(r["trajectory"])==1
            a=float(np.exp(grid[j+1]))
            assert abs(nxt["scale_factor"]/base["scale_factor"]/a-1.)<3e-12
            history.append(dict(block=j+1,relative_expansion=a,T_visible_GeV=nxt["T_visible_GeV"],
                                energy_ledger=r["energy_ledger"],count_ledger=r["count_ledger"],
                                accepted_microsteps=r["acceptances"],rejected_microsteps=r["rejections"],
                                exterior=r["exterior"],endpoint=r["endpoint"]))
            cur=nxt
            state=dict(plan=plan,plan_sha256=plan_hash,blocks_finished=j+1,blocks_total=len(grid)-1,
                       snapshot=cur,snapshot_sha256=hashlib.sha256(bytes_of(cur)).hexdigest(),ledger=history)
            atomic(saved,state)
            print(json.dumps(dict(blocks_finished=j+1,blocks_total=len(grid)-1,
                                  T_visible_GeV=cur["T_visible_GeV"],accepted=r["acceptances"],
                                  rejected=r["rejections"])),file=sys.stderr,flush=True)
        if step_end==len(grid)-1:
            atomic(directory/"final-shape.json",dict(snapshots=[cur]))
        return dict(plan_sha256=plan_hash,blocks_finished=step_end,blocks_total=len(grid)-1,
                    complete=step_end==len(grid)-1,final_T_visible_GeV=cur["T_visible_GeV"])


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for x in ("eos","start","directory"):
        p.add_argument(x,type=Path)
    p.add_argument("--stop",type=float,required=True)
    p.add_argument("--steps",type=int,required=True)
    p.add_argument("--local-prompt-jacobian",action="store_true")
    p.add_argument("--max-blocks",type=int)
    a=p.parse_args()
    print(json.dumps(main(a.eos,a.start,a.directory,a.stop,a.steps,a.local_prompt_jacobian,a.max_blocks),
                     indent=2,allow_nan=False))
