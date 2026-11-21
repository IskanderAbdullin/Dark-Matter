"""Deterministic microphysical 5D weak-form benchmark, without external BEST.

Run: python kinetics/verify_relativistic_weak_form.py --output path.json.
This is reproducible numerics, NOT evidence of 86-eV relic or full SM.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
from relativistic_weak_form_2to2 import (
    weak_energy_moments, rest_blocked_heat_independent)


def audit():
    def run(n,cutoff,Tchi=.8,Te=1.,G=.001):
        x=weak_energy_moments(
            M=.1,Tchi=Tchi,Te=Te,G=G,gchi=4,ge=4,
            n_p=n,n_k=n,n_mu=n,n_cm=n,n_phi=n,cutoff=cutoff)
        return {k:x[k] for k in (
            "Q_chi_GeV5","Q_e_GeV5","entropy_GeV4","number_strong_GeV4",
            "energy_strong_GeV5","weak_vs_strong_heat_relative_discrepancy",
            "max_mass_shell_defect","quadrature_points")}
    cases={
        "n12_cut24":run(12,24.),
        "n16_cut24":run(16,24.),
        "n24_cut32":run(24,32.),
        "n28_cut32":run(28,32.),
        "n24_cut15_negative_control":run(24,15.),
        "equilibrium_Tchi_eq_Te":run(12,24.,Tchi=1.),
        "hot_chi_cooling_Tchi_1p2":run(12,24.,Tchi=1.2),
    }
    a=cases["n24_cut32"];b=cases["n28_cut32"]
    neg=cases["n24_cut15_negative_control"]
    assert abs(b["Q_chi_GeV5"]/a["Q_chi_GeV5"]-1)<1e-4
    assert abs(b["weak_vs_strong_heat_relative_discrepancy"])<1e-5
    assert abs(b["number_strong_GeV4"])<1e-12
    assert neg["weak_vs_strong_heat_relative_discrepancy"]<-.005
    assert cases["equilibrium_Tchi_eq_Te"]["Q_chi_GeV5"]==0
    assert cases["hot_chi_cooling_Tchi_1p2"]["Q_chi_GeV5"]<0
    return {
        "status":"deterministic numerical reference (not an external solver)",
        "model":"tree vector contact chi e, massless e, isotropic FD species",
        "units":"M,T in GeV; G GeV^-2; Q GeV^5; entropy rate GeV^4",
        "benchmark":{"M_GeV":.1,"Tchi_GeV":.8,"Te_GeV":1.,
                     "G_GeVminus2":.001,"gchi":4,"ge":4,"mu":0.},
        "scope":"5D full-Pauli weak moment + sampled Cchi(p), not arbitrary f(p), gamma_el, BEST, or 86-eV freeze-out",
        "cutoff_warning":"finite pmax=cutoff*Tchi, kmax=cutoff*Te; direct strong number leakage and Q mismatch at low cutoff",
        "cases":cases,
        "independent_rest_Pauli_heat_M1_T0p1_G0p2":
            rest_blocked_heat_independent(M=1.,Te=.1,G=.2,n_k=24,n_cm=24),
        "no_gamma_el_or_relic_or_class":True,
    }


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="")
    args=p.parse_args()
    result=audit()
    if args.output:
        target=Path(args.output)
        target.parent.mkdir(parents=True,exist_ok=True)
        target.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))
