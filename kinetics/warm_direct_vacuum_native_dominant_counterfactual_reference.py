"""Prepare high-abundance shape tests, never infer them from the source yield."""
import argparse,json
from pathlib import Path


def prepare(high,out):
    out.mkdir(parents=True,exist_ok=True)
    base=(high/"compact-physical.ini").read_text()
    cases=[]
    for label,frac in (("ten_percent",.1),("half",.5)):
        q={}
        for l in base.splitlines():
            a,b=l.split("=",1);q[a.strip()]=b.strip()
        q["omega_cdm"]=format(.12*(1-frac),".17g")
        q["omega_ncdm"]=format(.12*frac,".17g")
        q["root"]=str(out.resolve()/(label+"_"))
        (out/(label+".ini")).write_text("".join(f"{k} = {v}\n" for k,v in q.items()))
        cases.append(dict(label=label,imposed_rest_fraction=frac))
    return dict(scope="Only a native counterfactual; retains the direct-source shape but imposes an unrelated, "
                      "much higher fraction, keeping the same present dark density. Neither production "
                      "nor statistical viability are inferred.",cases=cases)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("extended",type=Path);p.add_argument("output",type=Path)
    a=p.parse_args();print(json.dumps(prepare(a.extended,a.output),indent=2,allow_nan=False))
