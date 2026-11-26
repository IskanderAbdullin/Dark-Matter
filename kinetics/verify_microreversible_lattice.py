"""Reproduce direct microreversible 1D transfer + 3D energy full-Pauli heat.

Default runs χ/e h=1,.5,.25 GeV, 22GeV energy cutoff and independent
reference numerical 5D CM quadrature values measured at n=28.
No external BEST and no 86-eV relic.
Use --with-fine to include h=.125 (larger, optional).
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import numpy as np
from scipy.special import expit
from microreversible_chie_lattice import (
    invariant_kernel,microreversible_lattice,sources,evolve)


def audit(with_fine=False):
    ref={".9":1.9851492233308662e-8,"1.1":-4.432408224385169e-8}
    cases={}
    for h in ((1.,.5,.25,.125) if with_fine else (1.,.5,.25)):
        n=round(22/h)
        g=microreversible_lattice(dE=h,nchi=n,ne=n,nq=12)
        E=g["Echi"];k=g["Ee"];fe=expit(-k)
        data={"dE_GeV":h,"nchi":n,"ne":n,"quartets":len(g["W"])}
        for label,Tc in [(".9",.9),("1.1",1.1)]:
            result=sources(g,expit(-E/Tc),fe)
            Q=result["Qchi_GeV5"]
            data["Tchi_"+label+"_Q_GeV5"]=Q
            data["Tchi_"+label+"_relative_to_independent_5D"]=Q/ref[label]-1
            assert abs(result["energy_residual_GeV5"])<1e-17
            assert result["entropy_GeV4"]>0
        cases[str(h)]=data
    assert abs(cases["0.5"]["Tchi_.9_relative_to_independent_5D"])<5e-4
    assert abs(cases["0.5"]["Tchi_1.1_relative_to_independent_5D"])<5e-4
    probe=[]
    for E,k,Ep,kp in [(1.,2.,2.,1.),(3.,4.,4.,3.),(2.,4.,1.,5.),(.7,1.2,1.4,.5)]:
        f=float(invariant_kernel(.1,.001,E,k,np.array([Ep]),np.array([kp]))[0])
        r=float(invariant_kernel(.1,.001,Ep,kp,np.array([E]),np.array([k]))[0])
        assert math.isclose(f,r,rel_tol=1e-12)
        probe.append({"in":[E,k],"out":[Ep,kp],"rate_per_energy_cube_GeV":f,
                      "reverse_relative_difference":(r-f)/f})
    g=microreversible_lattice(dE=.5,nchi=14,ne=17,nq=12)
    E=g["Echi"];k=g["Ee"]
    F=expit(-E/1.2);R=expit(-k)
    initial=sources(g,F,R)["Qchi_GeV5"]
    time=evolve(g,F,R,1e3)
    assert abs(time["Qchi_integrated_GeV4"]/1e3/initial-1)<.003
    assert abs(time["energy_residual_GeV4"])<1e-12
    assert time["entropy_gain_GeV3"]>0
    return {
     "model":"electron-only massless e± broken-EW vector contact full Pauli χe",
     "units":"M,T,dE GeV; G GeV^-2; Q GeV5; rate kernel GeV",
     "parameters":{"Mchi":.1,"T_e":1.,"G":.001,"gchi":4,"ge":4,
                   "cutoff_energies_GeV":22,"nq":12},
     "independent_5D_Gauss_p_k_mu_CM_n28_reference_GeV5":ref,
     "energy_grid_scan":cases,
     "four_independent_time_reversal_energy_pairs":probe,
     "nonthermal_small_dt":{
       "initial_Qchi_GeV5":initial,"dt_GeVminus1":1e3,
       "delta_Echi_GeV4":time["Qchi_integrated_GeV4"],
       "energy_residual_GeV4":time["energy_residual_GeV4"],
       "chi_N_residual_GeV3":time["Nchi_residual"],
       "electron_N_residual_GeV3":time["Ne_residual"],
       "entropy_gain_GeV3":time["entropy_gain_GeV3"]},
     "microscopic_forward_reverse_paired_without_averaging":True,
     "all_bins_energy_collocation_approximation":True,
     "not_uv_gauge_completed_or_cosmological_Bose_relic":True,
     "class_ready":False
    }


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",default="")
    parser.add_argument("--with-fine",action="store_true")
    opt=parser.parse_args()
    result=audit(opt.with_fine)
    if opt.output:
        path=Path(opt.output)
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))
