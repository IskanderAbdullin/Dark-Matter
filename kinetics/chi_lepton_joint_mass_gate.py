"""Joint χ mass and lepton-mass 5D weak Pauli error screen.

Real-neutral-φ equal-charge model requires both χ species to scatter;
this only tests prescribed arbitrary fχ/fℓ, not their actual abundances
or coupled Yukawa/SM/FRW evolution. Do not infer physical Γ_el.
"""
from __future__ import annotations
import json
from kinetics.independent_chi_lepton_5d_quadrature import evaluate


def case(M,m,shape,level):
    kw={
      "coarse":dict(cutoff=12.,segments=12,nr=4,nmu=6,nc=12,nphi=12),
      "mid":dict(cutoff=12.,segments=12,nr=5,nmu=8,nc=16,nphi=16),
      "fine":dict(cutoff=16.,segments=16,nr=5,nmu=12,nc=32,nphi=24),
    }
    if level not in kw:raise ValueError("level")
    return evaluate(shape,M=M,m=m,mv=.1,g=.001,
       radial_scheme="segmented_p",angular_scheme="log_transfer",
       **kw[level])


def audit(full=False):
    out=[]
    for M in (.1,.02):
      for m in (0.,.5):
        for shape in ("smooth","oscillatory"):
          levels=("coarse","mid","fine") if full else ("coarse","mid")
          rows=[]
          for level in levels:
            x=case(M,m,shape,level)
            rows.append({"level":level,"Qchi":x["Qchi"],
                "entropy":x["entropy"],
                "max_energy_defect":x["max_energy_defect"],
                "negative_micro_entropy_events":x["negative_entropy_events"]})
          out.append({"Mchi":M,"mlepton":m,"shape":shape,"rows":rows})
    return {"schema_version":1,"cases":out,
       "scope":"finite-vector full on-shell Pauli weak χℓ Q/H for synthetic occupations, 2 masses each; NOT physical coupled χ1/χ2+SM/φ/FRW",
       "class_ready":False}


if __name__=="__main__":
    import argparse
    p=argparse.ArgumentParser()
    p.add_argument("--full",action="store_true")
    a=p.parse_args()
    print(json.dumps(audit(full=a.full),indent=2))
