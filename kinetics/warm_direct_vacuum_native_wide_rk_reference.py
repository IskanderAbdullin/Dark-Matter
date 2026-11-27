"""Extend the exact trace-centered RK comparison to higher native k."""
import argparse,json
from pathlib import Path

MODES=(.1,.3,1.,1.5,2.)
KINDS=("trace","physical","amp4","heldtrace","heldphysical")


def run(source,output):
    output.mkdir(parents=True,exist_ok=True)
    loads=[]
    for tier,nq,lmax in (("compact",70,20),("resolved",96,26)):
        for kind in KINDS:
            if tier=="compact" and kind.startswith("held"):continue
            label=tier+"-"+kind
            cfg={}
            for line in (source/(kind+".ini")).read_text().splitlines():
                k,v=line.split("=",1);cfg[k.strip()]=v.strip()
            cfg.update(dict(ncdm_N_momentum_bins=nq,l_max_ncdm=lmax))
            cfg["k_output_values"]=",".join(f"{v*.6736:.12g}" for v in MODES)
            cfg["P_k_max_h/Mpc"]="2.12"
            cfg["root"]=str(output.resolve()/(label+"_"))
            cfg["evolver"]=0
            cfg["tol_perturbations_integration"]=9e-10
            (output/(label+".ini")).write_text("".join(f"{k} = {v}\n" for k,v in cfg.items()))
            loads.append(dict(label=label,nq=nq,lmax=lmax,omega_nd=float(cfg["omega_ncdm"])))
    return dict(modes_h_per_Mpc=MODES,rows=loads,scope="Early-start exact RK, two explicit momentum/hierarchy tiers. "
                "Low-mode initialization already matched separately; no automatic fluid closure.")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("from_early",type=Path);p.add_argument("output",type=Path)
    a=p.parse_args();print(json.dumps(run(a.from_early,a.output),indent=2,allow_nan=False))
