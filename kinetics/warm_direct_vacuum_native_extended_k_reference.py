"""Prepare exact-mode early-start native RK checks through 10 h/Mpc."""
import argparse,json
from pathlib import Path

MODES=(.3,1.,2.,3.,5.,8.,10.)

def run(early,out):
    out.mkdir(parents=True,exist_ok=True)
    records=[]
    for tier,q,l,tol in (("dense",120,40,9e-10),("compact",78,28,1.2e-9)):
        for name in ("trace","physical","amp4","heldtrace","heldphysical"):
            if tier=="compact" and name in ("amp4","heldtrace","heldphysical"):continue
            cfg={}
            for line in (early/(name+".ini")).read_text().splitlines():
                k,v=line.split("=",1);cfg[k.strip()]=v.strip()
            cfg.update(dict(ncdm_N_momentum_bins=q,l_max_ncdm=l,tol_perturbations_integration=tol))
            cfg["k_per_decade_for_pk"]=8
            cfg["k_output_values"]=",".join(f"{x*.6736:.12g}" for x in MODES)
            cfg["P_k_max_h/Mpc"]="10.2"
            tag=tier+"-"+name
            cfg["root"]=str(out.resolve()/(tag+"_"))
            (out/(tag+".ini")).write_text("".join(f"{k} = {v}\n" for k,v in cfg.items()))
            records.append(dict(name=tag,quadrature=q,lmax=l,tolerance=tol,omega_ncdm=float(cfg["omega_ncdm"])))
    return dict(modes_h_per_Mpc=MODES,cases=records,scope="Already-produced scalar and a stable CDM proxy at fixed native background; two exact-mode RK tiers, not nonlinear clustering.")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("early",type=Path);p.add_argument("out",type=Path)
    a=p.parse_args();print(json.dumps(run(a.early,a.out),indent=2,allow_nan=False))
