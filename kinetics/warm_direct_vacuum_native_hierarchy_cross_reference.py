"""Cross the momentum and hierarchy controls at fixed exact-k RK settings."""
import argparse,json
from pathlib import Path


def run(wide,d):
    d.mkdir(parents=True,exist_ok=True)
    v=[]
    for type_,nq,lmax in (("q_only",96,20),("l_only",70,26)):
        for arm in ("trace","physical"):
            label=type_+"-"+arm
            q={}
            p=wide/("compact-"+arm+".ini")
            for line in p.read_text().splitlines():
                k,x=line.split("=",1);q[k.strip()]=x.strip()
            q["ncdm_N_momentum_bins"]=nq
            q["l_max_ncdm"]=lmax
            q["root"]=str(d.resolve()/(label+"_"))
            (d/(label+".ini")).write_text("".join(f"{k} = {x}\n" for k,x in q.items()))
            v.append(dict(label=label,nq=nq,lmax=lmax))
    return dict(scope="Only change one numerical axis from the compact exact RK tier at a time. "
                     "Both numerical axes must still be combined in the resolved corner.",rows=v)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("wide",type=Path);p.add_argument("cross",type=Path)
    a=p.parse_args();print(json.dumps(run(a.wide,a.cross),indent=2,allow_nan=False))
