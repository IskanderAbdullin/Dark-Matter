"""Extend only the radial q sampling after the hierarchy cross."""
import argparse
from pathlib import Path

def prepare(wide,out,points):
    out.mkdir(parents=True,exist_ok=True)
    for arm in ("trace","physical"):
        q={}
        for l in (wide/("compact-"+arm+".ini")).read_text().splitlines():
            k,v=l.split("=",1);q[k.strip()]=v.strip()
        q["ncdm_N_momentum_bins"]=points
        q["root"]=str(out.resolve()/(arm+"_"))
        (out/(arm+".ini")).write_text("".join(f"{k} = {v}\n" for k,v in q.items()))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("wide",type=Path);p.add_argument("out",type=Path)
    p.add_argument("--points",type=int,default=168)
    a=p.parse_args();prepare(a.wide,a.out,a.points)
