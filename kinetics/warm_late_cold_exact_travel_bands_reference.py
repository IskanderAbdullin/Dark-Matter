"""Exact cold-source recoil folded to actual birth-dependent travel cuts."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from archived_collisionless_streaming_reference import Background
from kloe_pion_band_portal_reference import table
from late_radiation_background_reference import TabulatedRadiation
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import P
from warm_late_cold_distance_source_reference import q_limit
from warm_late_cold_scalar_prompt_partition_reference import counts
from warm_late_cold_visible_exact_recoil_fold_reference import local
from warm_late_cold_visible_kloe_recoil_reference import rates
from warm_light_occupation_stage_audit_reference import digest,read


def run(eos,archive,kloe,distance,piondist):
    bg=Background(eos)
    S0=read(archive);S=S0["snapshot"]
    assert S0["snapshot_sha256"]==digest(S)
    d=read(distance);pion=read(piondist)
    gd,cs,err=table(kloe)
    Tr=.05
    rad=TabulatedRadiation(eos,Tr)
    T0=S["T_visible_GeV"]
    a0=bg.scale(T0)
    sv0=2*np.pi**2/45*T0**3*np.exp(bg.gs(np.log(T0)))
    records=[]
    for D in (.3,.5,1.,2.):
        q,_=q_limit(bg,a0,D)
        held=float(counts(S,np.array([q/a0]),24)[0]/sv0)
        rows=[]
        for no,nt in ((8,12),(12,20)):
            x,w=roots_legendre(nt)
            aa,bb=np.log([.012,Tr])
            dln=(bb-aa)/2
            vals=[]
            qs=[]
            for T in np.exp((aa+bb)/2+dln*x):
                a=bg.scale(T)
                q,_=q_limit(bg,a,D)
                pc=q/a
                assert pc<P
                qs.append(q/bg.T0)
                Z=local(T,pc,no)
                W=rates(T,pc,gd,cs,err,no)
                s=2*np.pi**2/45*T**3*np.exp(bg.gs(np.log(T)))
                _,rho=rad.state(a/bg.scale(Tr))
                dt=(3+bg.gs.derivative()(np.log(T)))/3
                vals.append(np.r_[Z[:2],W]*dt/(s*np.sqrt(rho/3)/MPL))
            y=dln*(w@np.array(vals))
            rows.append(dict(time_order=nt,kinetic_order=no,lepton_FD=float(y[0]),lepton_MB=float(y[1]),
                             pion_MB=float(y[2]),pion_BE=float(y[3]),pion_linear_marginal=float(y[4]),
                             combined_central_over_held=float((y[0]+y[3])/held),
                             qcut_over_T0_range=[min(qs),max(qs)]))
        assert abs(rows[0]["combined_central_over_held"]/rows[-1]["combined_central_over_held"]-1)<.01
        if D==.3:
            r=d["rows"][2]
            pi=pion["cases"][2]
            assert r["max_path_Mpc"]==pi["max_distance_Mpc"]==D
            assert rows[-1]["combined_central_over_held"]<pi["combined_conditional_upper_over_held_at_high_R_1"]
        records.append(dict(max_actual_path_Mpc=D,held_distance_selected_yield=held,controls=rows))
    return dict(scope="Numerical 12–50 MeV Fermi-Dirac charged leptons plus measured charged pion vacuum window, "
                 "using the actual injection-time-dependent reference path cutoff. No dark scattering or complete hadronic source.",
                archive=S0["snapshot_sha256"],KLOE_sha256=hashlib.sha256(kloe.read_bytes()).hexdigest(),rows=records)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","archive","kloe","distance","piondist"):
        p.add_argument(k,type=Path)
    a=p.parse_args();print(json.dumps(run(a.eos,a.archive,a.kloe,a.distance,a.piondist),indent=2,allow_nan=False))
