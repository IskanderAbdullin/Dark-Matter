"""Reproducible compact 5D χe equilibrium bin-projection bias experiment.

No huge matrices are saved. All values refer to electron-only vector EFT,
massless e±, equilibrium χ and e at T=1 GeV. Not a full kinetic solver.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.special import expit
from microphysical_equilibrium_graph import (
    graph,graph_instantaneous_heat,calibrated_pair_rates,local_pair_implicit_step)

FIELDS=("D_true_all","D_true_inside","D_graph",
        "jensen_relative","unpaired_asymmetry",
        "raw_stationary_relative_L1",
        "raw_FD_heat_without_reversible_closure_GeV5",
        "outgoing_out_of_support_weight_fraction")


def audit():
    configs={
        "chi16_small":dict(nchi=16,nk=12,nmu=10,ncm=10,nphi=10),
        "chi24_small":dict(nchi=24,nk=12,nmu=10,ncm=10,nphi=10),
        "chi64_angular12":dict(nchi=64,nk=12,nmu=12,ncm=12,nphi=12),
        "chi64_angular16":dict(nchi=64,nk=16,nmu=16,ncm=16,nphi=16),
        "chi64_angular24":dict(nchi=64,nk=24,nmu=24,ncm=24,nphi=24),
        "chi128_nk20_ang16":dict(nchi=128,nk=20,nmu=16,ncm=16,nphi=16)
    }
    result={}
    for key,config in configs.items():
        g=graph(**config)
        result[key]={k:float(g[k]) for k in FIELDS}
        result[key]["config"]=config
        if key=="chi64_angular16":
            E=g["E"];f=expit(-E)
            Q,rate=graph_instantaneous_heat(g,f)
            eps=1e-4
            Qp,_=graph_instantaneous_heat(g,expit(-E+eps*E))
            Qm,_=graph_instantaneous_heat(g,expit(-E-eps*E))
            result[key]["closed_FD_equilibrium_Q_GeV5"]=float(Q)
            result[key]["closed_equilibrium_net_max_GeV4"]=float(np.max(np.abs(rate)))
            result[key]["centered_slope_GeV4"]=float((Qp-Qm)/(2*eps))
            result[key]["relative_to_minus_projected_D"]=float(
                (Qp-Qm)/(2*eps)/(-g["D_graph"])-1.)
            i,j,_=max(calibrated_pair_rates(g),key=lambda x:x[2])
            fi=float(expit(-E[i]+.02*E[i]));fj=float(expit(-E[j]+.02*E[j]))
            step=local_pair_implicit_step(g,i,j,fi,fj,10.,10.,1e6)
            result[key]["frozen_reference_T_finite_bath_pair"]={
                "i":int(i),"j":int(j),
                "Gamma_GeV4":float(step["Gamma_at_reference_T"]),
                "number_residual":float(step["particle_residual"]),
                "energy_residual":float(step["energy_residual"]),
                "entropy_gain":float(step["entropy_delta"]),
                "delta_T_GeV":float(step["T"]-1.),
                "is_actual_temperature_dependent_electron_integral":False
            }
    ref=3.027899802221842e-7
    assert abs(result["chi128_nk20_ang16"]["D_true_all"]/ref-1)<.001
    assert result["chi16_small"]["jensen_relative"]>.07
    assert result["chi64_angular24"]["unpaired_asymmetry"]<.008
    assert abs(result["chi64_angular16"]["relative_to_minus_projected_D"])<1e-6
    return {
       "source":"5D Lorentz-invariant physical χe vector contact matrix element, full outgoing electron and χ Pauli, equilibrium FD",
       "benchmark":{"Mchi_GeV":.1,"T_GeV":1.,"G_GeVminus2":.001,
                    "gchi":4,"ge":4,"pmax_GeV":22,"kmax_GeV":22},
       "independent_full_5D_thermal_D11_GeV4":ref,
       "scopes":"Raw projected conductances are physically integrated but nonreversible on finite bins; reversible closure is a declared APPROXIMATION at reference T, not exact nonlinear scattering.",
       "results":result,
       "class_ready":False,"physical_86eV_relic":False
    }


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="")
    args=p.parse_args()
    result=audit()
    if args.output:
        path=Path(args.output)
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2,ensure_ascii=False))
