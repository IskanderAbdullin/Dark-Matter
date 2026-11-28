"""Exact durable restart state for a locked comoving evolution."""
from contextlib import contextmanager
import fcntl
from pathlib import Path

import numpy as np

from frw_exterior_cohort_reference import Exterior
from warm_coupled_serial_restart_reference import atomic,sha
from warm_light_occupation_stage_audit_reference import digest,read


SOURCE_FILES=(
    "warm_finite_parent_comoving_expansion_reference.py",
    "warm_comoving_durable_restart_reference.py",
    "warm_finite_comoving_sdirk_reference.py",
    "warm_comoving_parent_chart_reference.py",
    "warm_comoving_parent_competing_reference.py",
    "warm_late_finite_parent_prompt_reference.py",
    "warm_compact_absolute_prompt_reference.py",
    "joint_prompt_two_moment_evolution_reference.py",
    "joint_resonant_centered_stiff_reference.py",
    "joint_resonant_frw_split_reference.py",
    "frw_exterior_cohort_reference.py",
    "warm_pair_common_fermion_chart_reference.py",
    "late_radiation_background_reference.py",
)


def key(eos,paths,options,qhash,whash):
    directory=Path(__file__).resolve().parent
    return digest(dict(version=1,eos=sha(eos),input_tapes=[sha(p) for p in paths],
                       numeric_code={p:sha(directory/p) for p in SOURCE_FILES},
                       options=options,parent_q=qhash,parent_weights=whash))


@contextmanager
def guard(folder):
    folder.mkdir(exist_ok=True,parents=True)
    with (folder/".restart-lock").open("a") as f:
        fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB)
        yield


def seal(filename,lineage,last,hv,work,ex,outside,exchange,max_pair_defect,row_digests):
    payload=dict(schema="locked-parent-frw-v1",key=lineage,accepted_blocks=last,
                 parent_occupations=np.asarray(hv).tolist(),pressure_work=float(work),
                 exterior=dict(a=float(ex.a),trace=float(ex.trace),injected=float(ex.injected),
                               cohorts=ex.cohorts),
                 outside=outside.tolist(),exchange=exchange.tolist(),
                 max_pair_defect=float(max_pair_defect),accepted_row_digests=row_digests)
    atomic(filename,dict(payload=payload,sha256=digest(payload)))


def restore(filename,lineage,scales,initial_hash,Q,H0):
    s=read(filename)
    x=s["payload"]
    assert s["sha256"]==digest(x) and x["schema"]=="locked-parent-frw-v1"
    assert x["key"]==lineage
    k=x["accepted_blocks"]
    assert isinstance(k,int) and 1<=k<len(scales)
    v=np.asarray(x["parent_occupations"])
    assert v.shape==Q.shape and np.all(np.isfinite(v)) and np.all((0.<=v)&(v<=1.))
    ex=Exterior()
    for key_ in ("a","trace","injected","cohorts"):
        setattr(ex,key_,x["exterior"][key_])
    assert ex.a==scales[k]
    assert np.all(np.isfinite(ex.array()))
    assert all(1.<=item[0]<=ex.a and item[5]>0. for item in ex.cohorts)
    ex.result()
    rows=[]
    ld=x["accepted_row_digests"]
    assert len(ld)==k
    for j,ph in enumerate(ld,1):
        y=read(filename.parent/f"step-{j}.json")
        r=y["row"]
        assert ph==digest(y)
        assert r["daughter_hash"]==digest(y["daughter_snapshot"])
        assert r["c"]==scales[j] and r["initial_daughter_hash"]==initial_hash
        rows.append(r)
    assert H0(v)==rows[-1]["parent_hash"] or rows[-1]["parent_hash"] is None
    y=read(filename.parent/f"step-{k}.json")
    return x,rows,y["daughter_snapshot"],ex,v
