"""Compact measured physical χe 5D -> Pauli quartet convergence control.

Reproduce: python kinetics/verify_paired_energy_quartets.py --output result.json
This is a finite energy lattice with explicit pair-averaged raw scattering
rates, not an unmodified continuum Boltzmann operator or a DM relic.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.special import expit
from paired_energy_quartets import quartets,sources,implicit_one_quartet

def audit():
    configs={
      "chi24_e48_de0p5_ang10":dict(nchi=24,ne=48,dE=.5,nmu=10,ncm=10,nphi=10),
      "chi36_e54_de0p5_ang10":dict(nchi=36,ne=54,dE=.5,nmu=10,ncm=10,nphi=10),
      "chi30_e30_de0p4_ang12":dict(nchi=30,ne=30,dE=.4,nmu=12,ncm=12,nphi=12),
      "chi30_e30_de0p4_ang16":dict(nchi=30,ne=30,dE=.4,nmu=16,ncm=16,nphi=16),
    }
    result={}
    for name,params in configs.items():
        g=quartets(**params)
        E=g["Echi"];k=g["Ee"];electron=expit(-k)
        cold=sources(g,expit(-E/.9),electron)
        hot=sources(g,expit(-E/1.1),electron)
        equilibrium=sources(g)
        raw=sources(g,raw=True)
        result[name]={
          "parameters":params,
          "D_projected_quartet_GeV4":g["D_quartet_discrete"],
          "D_independent_5D_event_sample_GeV4":g["D_continuous_5D"],
          "raw_event_Pauli_weighted_reverse_asymmetry":g["W_eq_Pauli_asym_L1"],
          "FD_weighted_out_of_support_fraction":g["omitted_FD_weight_fraction"],
          "raw_unweighted_out_of_support_fraction":g["omitted_raw_weight_fraction"],
          "raw_FD_spurious_heat_GeV5":raw["Qchi"],
          "paired_FD_heat_GeV5":equilibrium["Qchi"],
          "cold_chi_0p9_GeV5":cold["Qchi"],
          "hot_chi_1p1_GeV5":hot["Qchi"],
          "hot_entropy_GeV3_per_time_GeV4":hot["entropy"],
          "total_energy_exchange_residual":hot["energy_total"],
        }
        assert abs(equilibrium["Qchi"])<1e-18
        assert abs(hot["Qchi"]+hot["Qe"])<1e-18
        assert hot["entropy"]>0 and cold["entropy"]>0
    g=quartets(nchi=30,ne=30,dE=.4,nmu=8,ncm=8,nphi=8)
    E=g["Echi"];k=g["Ee"];W=g["W"];f=expit(-E);e=expit(-k)
    I,A,J,B=np.unravel_index(
        np.argmax(W*(E[:,None,None,None]<E[None,None,:,None])),W.shape)
    step=implicit_one_quartet(g,int(I),int(A),int(J),int(B),
        float(.7*f[I]),float(1.1*e[A]),float(1.1*f[J]),float(.9*e[B]),1e8)
    assert abs(step["energy_residual_GeV4"])<1e-12
    assert step["entropy_change_GeV3"]>=0
    benchmark_cold=1.9852809333031555e-8
    benchmark_hot=-4.432526721786192e-8
    x=result["chi36_e54_de0p5_ang10"]
    assert abs(x["cold_chi_0p9_GeV5"]/benchmark_cold-1)<.025
    assert abs(x["hot_chi_1p1_GeV5"]/benchmark_hot-1)<.025
    return {
        "status":"physical 5D χe matrix element / full Pauli, with explicit paired finite-lattice rate approximation and a verified implicit quartet substep",
        "model":{"M_GeV":.1,"Te_GeV":1.,"G_GeVminus2":.001,
                 "gchi":4,"ge":4,"massless_e":True},
        "independent_continuum_5D_weak_reference":{
          "Tchi_0p9_Q_GeV5":benchmark_cold,
          "Tchi_1p1_Q_GeV5":benchmark_hot,
          "method":"separate Gauss 5D p/k/CM numerical source on [0,22]GeV, n14, no bins; not an external BEST execution"
        },
        "experiments":result,
        "implicit_stiff_pair":{"quartet_indices":[int(I),int(A),int(J),int(B)],
         "rate_GeV4":step["pair_averaged_5D_W_GeV4"],
         "transfer_number_GeV3":step["transfer_number_GeV3"],
         "energy_residual_GeV4":step["energy_residual_GeV4"],
         "chi_charge_residual":step["chi_charge_residual"],
         "electron_charge_residual":step["electron_charge_residual"],
         "entropy_change_GeV3":step["entropy_change_GeV3"],
         "is_exact_continuum_kernel":False},
        "negative_result":"Raw physical event bin tensor W_ia,jb is not reversible; pair averaging restores equilibrium but modifies finite-bin dynamics, FD-weighted truncation nonzero",
        "not_full_global_time_evolution":True,
        "gamma_el_or_86eV_relic_validated":False,
        "class_ready":False
    }

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",default="")
    args=parser.parse_args()
    result=audit()
    if args.output:
        path=Path(args.output)
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2,ensure_ascii=False))
