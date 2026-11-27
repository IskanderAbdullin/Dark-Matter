"""Generate reproducible quadrature/hierarchy tiers for the thermostatic exact-mode comparison."""
import argparse
import json
from pathlib import Path

TIERS={
    "linear78_cap45":(3,78,45,28,"1.2e-09"),
    "linear136_cap45":(3,136,45,40,"7e-10"),
    "linear136_cap18":(3,136,18,40,"7e-10"),
    "transform100":(2,100,45,40,"7e-10"),
    "transform160":(2,160,45,48,"4e-10"),
}
NAMES=("direct50_trace","direct50_absolute","direct100_trace","direct100_absolute","direct100_matched50")


def run(original,out,tier):
    strategy,n,cap,lmax,tol=TIERS[tier]
    out.mkdir(parents=True,exist_ok=True)
    rows=[]
    for label in NAMES:
        c=dict((a.strip(),b.strip()) for a,b in (v.split("=",1) for v in (original/(label+".ini")).read_text().splitlines()))
        c.update(ncdm_quadrature_strategy=str(strategy),ncdm_N_momentum_bins=str(n),
                 ncdm_maximum_q=str(cap),l_max_ncdm=str(lmax),tol_perturbations_integration=tol,
                 root=str((out/(label+"_")).resolve()))
        (out/(label+".ini")).write_text("".join(f"{k} = {v}\n" for k,v in c.items()))
        rows.append(dict(label=label,omega_ncdm=float(c["omega_ncdm"]),gauge=c["gauge"]))
    return dict(tier=tier,strategy=strategy,n=n,qcap=cap if strategy==3 else None,lmax=lmax,
                integration_tolerance=tol,cases=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("original",type=Path)
    p.add_argument("out",type=Path)
    p.add_argument("tier",choices=tuple(TIERS))
    a=p.parse_args()
    print(json.dumps(run(a.original,a.out,a.tier),indent=2,allow_nan=False))
