"""Reproducible reduced full χ/e paired-quartet kinetic history with cutoff warning.

Run: python kinetics/verify_paired_quartet_network.py --output out.json
No physical full Boltzmann/SM/Bose cosmological prediction.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
from scipy.special import expit
from paired_energy_quartets import quartets,sources
from paired_quartet_network import evolve


def audit():
    g=quartets(nchi=14,ne=17,dE=.5,nmu=8,ncm=8,nphi=8)
    f=expit(-g["Echi"]/1.2);e=expit(-g["Ee"])
    Q=sources(g,f,e)["Qchi"]
    cases={}
    for n in (1,4,16,64):
        forward=evolve(g,f,e,2e6,substeps=n)
        backward=evolve(g,f,e,2e6,substeps=n,reverse=True)
        cases[str(n)]={
           "delta_Echi_GeV4":forward["delta_Echi"],
           "delta_Ee_GeV4":forward["delta_Ee"],
           "total_energy_residual_GeV4":forward["total_energy_residual"],
           "Nchi_residual_GeV3":forward["Nchi_residual"],
           "Ne_residual_GeV3":forward["Ne_residual"],
           "total_Fermi_entropy_gain_GeV3":forward["entropy_gain"],
           "forward_reverse_order_heat_discrepancy_GeV4":
               abs(forward["delta_Echi"]-backward["delta_Echi"]),
           "number_of_physical_rate_quartets":forward["edges"]
        }
        assert forward["entropy_gain"]>0
        assert abs(forward["total_energy_residual"])<1e-12
        assert abs(forward["Nchi_residual"])<1e-12
        assert abs(forward["Ne_residual"])<1e-12
    tiny={}
    for dt in (1e3,1e4,1e5):
        x=evolve(g,f,e,dt,substeps=1)
        tiny[str(dt)]={
         "finite_time_dEchi_dt_GeV5":x["delta_Echi"]/dt,
         "relative_error_to_independent_instantaneous_Q":
           x["delta_Echi"]/dt/Q-1
        }
    assert (cases["16"]["forward_reverse_order_heat_discrepancy_GeV4"] <
            cases["4"]["forward_reverse_order_heat_discrepancy_GeV4"])
    common=evolve(g,expit(.15-g["Echi"]),expit(-.2-g["Ee"]),1e7,substeps=4)
    assert abs(common["delta_Echi"])<1e-12
    return {
      "status":"reproducible all-quartet finite-lattice Pauli χ/e dynamics with PAIR-AVERAGED physical 5D rates, NOT unmodified continuum microphysics",
      "benchmark":{"M_GeV":.1,"T_e_initial_GeV":1.,"Tchi_initial_GeV":1.2,
                   "G_GeVminus2":.001,"gchi":4,"ge":4,
                   "nchi":14,"ne":17,"deltaE_GeV":.5,
                   "angular_n":8,"duration_GeVminus1":2e6},
      "initial_instantaneous_Qchi_GeV5":Q,
      "time_step_and_order_scan":cases,
      "small_dt_reference":tiny,
      "common_T_stationarity_with_independent_mu":True,
      "limitations":[
          "All W are explicitly PAIR-AVERAGED finite-bin projections of physical 5D chi-e scattering; raw W is not reversible",
          "No electron mass, neutrino, UV mediator or SM gstar, no FRW or actual Bose-stimulated 1<->2",
          "A positive entropy, charges and energy do not imply dt/angle/grid convergence",
          "No 86-eV DM relic, BEC, CLASS, or validated kinetic gamma"
      ],
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
        path.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))
