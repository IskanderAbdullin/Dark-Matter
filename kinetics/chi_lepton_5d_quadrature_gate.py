"""Independently vary one quadrature axis of the nonthermal χℓ 5D weak rate.

Positive microscopic entropy and exact Qχ=-Qℓ do NOT establish finite
Gaussian-grid convergence. These prescribed synthetic distributions
are NOT early-universe kinetic solutions or physical Γ_el.
"""
from __future__ import annotations
import argparse
import json
from kinetics.independent_chi_lepton_5d_quadrature import evaluate


def audit(full=False):
    scans={}
    basis={"nr":8,"nmu":8,"nc":24,"nphi":16,"cutoff":12.}
    variations={
      "incoming_radial":("nr",[4,6,8,10,12,16] if full else [4,6,8,10]),
      "incoming_angle":("nmu",[4,6,8,10,12,16] if full else [4,6,8,10]),
      "outgoing_cm_polar":("nc",[8,12,16,24,32,48] if full else [8,12,16,24,32]),
      "outgoing_azimuth":("nphi",[8,12,16,24,32] if full else [8,12,16,24]),
      "kinetic_cutoff":("cutoff",[8.,12.,16.,20.] if full else [8.,12.,16.]),
    }
    for shape in ("smooth","oscillatory"):
      scans[shape]={}
      for axis,(key,values) in variations.items():
        scans[shape][axis]=[]
        for v in values:
          kw=dict(basis)
          kw[key]=v
          sample=evaluate(shape,**kw)
          scans[shape][axis].append({key:v,"Qchi":sample["Qchi"],
            "entropy":sample["entropy"],
            "negative_micro_entropy_events":sample["negative_entropy_events"],
            "energy_defect":sample["max_energy_defect"]})
    return {"schema_version":1,
       "scope":"synthetic incoming-integrated nonthermal Pauli 2to2 finite-vector weak heat/entropy; NOT continuum- or cosmology-certified",
       "fixed":{"Mchi":.1,"mlepton":0,"g_effective":.001,"mV":.1,
                "gspin_chi":4,"gspin_lepton":4},
       "reference_parameters":basis,"axes":scans,
       "physical_dm_relic_validated":False}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--full",action="store_true")
    args=parser.parse_args()
    print(json.dumps(audit(full=args.full),indent=2))
