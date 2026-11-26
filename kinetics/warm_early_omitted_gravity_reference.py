"""Restore only an instantaneous bound on omitted old exterior gravity."""
import argparse
import json
from pathlib import Path

import numpy as np

from warm_archived_expansion_exterior_envelope_reference import MASSES, from_ledger, history, retained
from warm_early_support_time_interaction_reference import read, sha
from frw_exterior_restart_reference import import_restart
from late_radiation_background_reference import TabulatedRadiation


def run(eos,budget,production,seeds,q4,q8):
    A=read(budget)["cases"]
    H=sha(eos)
    ans=[]
    for j,index in enumerate((16,32)):
        s=read(seeds[j])["snapshots"][-1]
        p=read(production[j])["rows"][-1]
        C=[read(x[j]) for x in (q4,q8)]
        for n,r in zip((4,8),C,strict=True):
            assert r["blocks_finished"]==r["blocks_total"]==n
            assert r["plan"]["eos_sha256"]==H and r["plan"]["source_sha256"]==sha(seeds[j])
        a=s["scale_factor"]
        assert abs(s["T_visible_GeV"]/p["final_visible_temperature_GeV"]-1.)<1e-12
        norm=np.asarray(A[j]["retained_fixed_Omega_h2"])/(MASSES*np.asarray(s["represented_species_count_a3"]))
        assert abs(norm[0]/norm[1]-1.)<3e-12
        N=np.asarray(A[j]["exterior_fixed_Omega_h2"])/(MASSES*norm)
        old=[from_ledger(p["visible_bath_feedback"],a**4,N,a)]
        assert abs(C[0]["snapshot"]["T_visible_GeV"]/C[1]["snapshot"]["T_visible_GeV"]-1.)<1e-12
        data=[]
        for n,r in zip((0,4,8),(s,C[0]["snapshot"],C[1]["snapshot"]),strict=True):
            rho=retained(r,44)[:,1].sum()
            rr=retained(r,72)[:,1].sum()
            assert abs(rho/rr-1.)<4e-10
            T=r["T_visible_GeV"]
            rad=TabulatedRadiation(eos,T).state(1.)[1]
            native_base=rad+rr
            if n:
                new=import_restart(r["positive_exterior_restart"])
                native_base+=new.energy(1.)
            ext,pressure=history(old,r["scale_factor"])
            d=ext/native_base
            dr=d/(np.sqrt(1+d)+1.)
            data.append(dict(outer_steps=n,photon_T_GeV=T,approximate_retained_and_new_exterior_over_radiation=native_base/rad-1.,
                             old_ext_rho_interval_GeV4=ext.tolist(),old_ext_pressure_interval_GeV4=pressure.tolist(),
                             relative_old_ext_added_H_interval=dr.tolist(),old_ext_added_energy_over_radiation=(ext/rad).tolist()))
        ans.append(dict(initial_support=index,normalized_old_exterior_inclusive_rest_omega=A[j]["exterior_fixed_Omega_h2"],data=data))
    return dict(cases=ans,EOS_sha256=H,input_record_sha256=[sha(p) for p in (budget,*production,*seeds,*q4,*q8)],
        scope="Instantaneous bracket on the historical exterior gravity excluded from the short restarted model, assuming those older tails freely redshift and never collide or return. Uses their known archive number and energy bounds, not a guessed momentum distribution. New exterior gravity was included in the numerical step. This is not a trajectory rerun or a bound on the induced reaction or abundance error.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos",type=Path)
    p.add_argument("budget",type=Path)
    for k in ("production","seeds","q4","q8"):
        p.add_argument("--"+k,type=Path,nargs=2,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.budget,a.production,a.seeds,a.q4,a.q8),indent=2,allow_nan=False))
