"""Use an explicit momentum rescaling for an illustrative uniform entropy dilution."""
import argparse,json
from pathlib import Path
from warm_light_occupation_stage_audit_reference import read


def run(high,abundance,out):
    A=read(abundance)
    D=A["dilution_factoring_old_fraction_to_reference_total"]
    eta=D**(-1/3)
    assert 0.<eta<1.
    out.mkdir(parents=True,exist_ok=True)
    rows=[]
    for label in ("trace","physical","amp4"):
        arm="dense" if label=="amp4" else "compact"
        s={}
        for l in (high/(arm+"-"+label+".ini")).read_text().splitlines():
            x,v=l.split("=",1);s[x.strip()]=v.strip()
        for k,v in dict(ncdm_N_momentum_bins=78,l_max_ncdm=28,tol_perturbations_integration="1.2e-9").items():
            s[k]=v
        s["T_ncdm"]=format(eta,".17g")
        s["root"]=str(out.resolve()/(label+"_"))
        (out/(label+".ini")).write_text("".join(f"{k} = {v}\n" for k,v in s.items()))
        rows.append(dict(label=label,omega_scalar=float(s["omega_ncdm"])))
    return dict(scope="Only a shape counterfactual: changing the final thermal momentum unit by D^(-1/3) "
                     "as an illustrative homogeneous visible entropy addition. It does not model the "
                     "entropy-release dynamics, radiation sector, perturbations or physical feasibility.",
                dilution=D,momentum_scaling=eta,cases=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("higher","normalization","out"):
        p.add_argument(k,type=Path)
    a=p.parse_args();print(json.dumps(run(a.higher,a.normalization,a.out),indent=2,allow_nan=False))
