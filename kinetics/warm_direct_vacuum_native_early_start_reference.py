"""Match initial times by forcing an earlier common opacity threshold."""
import argparse,json
from pathlib import Path

def run(source,dest,limit,templates=None):
    dest.mkdir(parents=True,exist_ok=True)
    labels=["base","trace","physical","amp4"]
    for name in labels:
        q={}
        for line in (source/(name+".ini")).read_text().splitlines():
            k,v=line.split("=",1);q[k.strip()]=v.strip()
        q["start_small_k_at_tau_c_over_tau_h"]=format(limit,".12g")
        q["root"]=str(dest.resolve()/(name+"_"))
        (dest/(name+".ini")).write_text("".join(f"{k} = {v}\n" for k,v in q.items()))
    if templates is not None:
        for tag in ("trace","physical"):
            x=(dest/(tag+".ini")).read_text().replace("direct-sig07.dat","held-sig07.dat")
            x=x.replace(str(dest.resolve()/(tag+"_")),str(dest.resolve()/("held"+tag+"_")))
            (dest/("held"+tag+".ini")).write_text(x)
        genesis=(dest/"amp4.ini").read_text()
        for n in ("direct3","p12","thermal3"):
            lines=[]
            for line in genesis.splitlines():
                if line.startswith("omega_cdm"):line="omega_cdm = 0.11988"
                if line.startswith("omega_ncdm"):line="omega_ncdm = 0.00012"
                if line.startswith("root"):line=f"root = {dest.resolve()/n}_"
                if n=="p12" and line.startswith("ncdm_psd_filenames"):
                    line=f"ncdm_psd_filenames = {(templates/'fd-p12-mix.dat').resolve()}"
                if n=="thermal3":
                    if line.startswith(("ncdm_psd_filenames","use_ncdm_psd_files")):continue
                    if line.startswith("T_ncdm"):
                        line=next(s for s in (templates/"std-thermal3.ini").read_text().splitlines() if s.startswith("T_ncdm"))
                lines.append(line)
            (dest/(n+".ini")).write_text("\n".join(lines)+"\n")
    spec=json.loads((source/"spec.json").read_text())
    spec["common_opacity_start_threshold"]=limit
    (dest/"spec.json").write_text(json.dumps(spec,indent=2)+"\n")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("source",type=Path);p.add_argument("dest",type=Path)
    p.add_argument("--opacity",type=float,default=3e-6)
    p.add_argument("--templates",type=Path)
    a=p.parse_args();run(a.source,a.dest,a.opacity,a.templates)
