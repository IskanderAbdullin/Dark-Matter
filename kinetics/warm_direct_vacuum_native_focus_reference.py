"""Prepare a narrowed independent RK and exact-probe native comparison."""
import argparse,json
from pathlib import Path

K=(.002,.01,.03,.1,.3,1.)


def prepare(input_dir,output):
    output.mkdir(parents=True,exist_ok=True)
    result=[]
    for label in ("base","trace","physical","amp4"):
        s={}
        for line in (input_dir/("std-"+label+".ini")).read_text().splitlines():
            k,v=line.split("=",1)
            s[k.strip()]=v.strip()
        s.update({"P_k_max_h/Mpc":"1.15","k_output_values":",".join(f"{q*.6736:.12g}" for q in K),
                  "k_per_decade_for_pk":"18","evolver":"0","tol_perturbations_integration":"9e-10",
                  "root":str(output.resolve()/(label+"_"))})
        if "N_ncdm" in s:
            s.update({"ncdm_N_momentum_bins":"96","l_max_ncdm":"24",
                      "ncdm_psd_filenames":str((input_dir/"direct-sig07.dat").resolve())})
        (output/(label+".ini")).write_text("".join(f"{k} = {v}\n" for k,v in s.items()))
        result.append(dict(tag=label,omega_cdm=float(s["omega_cdm"]),omega_ncdm=float(s.get("omega_ncdm","0"))))
    return dict(scope="Independent narrowed native RK response. Output probes lie inside 1.15 h/Mpc support; "
                      "do not label this as the original higher-tier 128q/l32 ladder.",
                k_requested_h_per_Mpc=K,tolerance=9e-10,quadrature_nodes=96,max_multipole=24,rows=result)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("input",type=Path);p.add_argument("output",type=Path)
    a=p.parse_args();print(json.dumps(prepare(a.input,a.output),indent=2,allow_nan=False))
