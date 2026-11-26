"""Sequential matched-block calibration of the two warm prompt derivatives."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

from warm_exterior_restart_audit_reference import shell


def load(p):
    return json.loads(p.read_text())


_CHILD = """
import json,resource,subprocess,sys,time
from pathlib import Path
out=Path(sys.argv[1])
t=time.perf_counter()
r=subprocess.run(sys.argv[2:])
u=resource.getrusage(resource.RUSAGE_CHILDREN)
out.write_text(json.dumps(dict(exit_code=r.returncode,wall_seconds=time.perf_counter()-t,
    user_seconds=u.ru_utime,kernel_seconds=u.ru_stime,peak_RSS_KiB=u.ru_maxrss))+"\\n")
raise SystemExit(r.returncode)
"""


def measure(argv,out,log):
    note=out.with_suffix(".resources")
    launcher=[sys.executable,"-W","error","-c",_CHILD,str(note),*argv]
    with out.open("w") as f, log.open("w") as stderr:
        subprocess.run(launcher,stdout=f,stderr=stderr,check=True)
    return load(note)

def run(eos,shape,out,expansions):
    out.mkdir(parents=True,exist_ok=True)
    rows=[]
    exe=Path(__file__).resolve().parent/"coupled_warm_pair_continuation_reference.py"
    for b in expansions:
        mm=[]
        for method in ("standard","local"):
            prefix=out/f"{b:g}-{method}"
            shape_out=Path(str(prefix)+"-shape.json")
            record=Path(str(prefix)+".json")
            cmd=[sys.executable,"-W","error",str(exe),str(eos),str(shape),"--stop",str(b),
                 "--steps","1","--save-snapshots",str(shape_out),"--progress-every","1"]
            if method=="local":
                cmd.append("--local-prompt-jacobian")
            rm=measure(cmd,record,Path(str(prefix)+".err"))
            mm.append((shell(load(shape_out)["snapshots"][-1]),load(record)["rows"][-1],rm))
        A,B=mm[0][0],mm[1][0]
        err=float(np.max(np.sum(abs(A-B),axis=2)/np.sum(abs(A),axis=2)))
        assert err<1e-9
        rows.append(dict(relative_expansion=b,shell_moment_L1_relative=err,
                         accepted_rejected=[[m[1]["acceptances"],m[1]["rejections"]] for m in mm],
                         resource_observations=[m[2] for m in mm]))
    return dict(rows=rows,scope="Same short blocks, not a full warm convergence or hardware guarantee.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for k in ("eos","shape","out"):
        p.add_argument(k,type=Path)
    p.add_argument("--expansions",type=float,nargs="+",default=[1.0002,1.03])
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.shape,a.out,a.expansions),indent=2,allow_nan=False))
