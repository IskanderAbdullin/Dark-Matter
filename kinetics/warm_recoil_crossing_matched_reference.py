"""Follow the cold-light scalar recoil gate on matched q32 and q64 trajectories."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from joint_prompt_two_moment_evolution_reference import Collision
from late_radiation_background_reference import TabulatedRadiation
from three_species_decay_onsager_reference import mp,width
from warm_light_occupation_stage_audit_reference import digest,read
from warm_near_rest_light_inverse_decay_reference import E0,freeze,window
from warm_cold_prompt_support_grid_reference import independent_first_average
from warm_prompt_reentry_shape_audit_reference import get_mom


def run(eos,D,work,limit):
    A=read(D/"warm-history-lower-q32-16-terminal-2026-11-26.json.gz")[-1]
    S0=A["snapshot"]
    assert digest(S0)==A["snapshot_sha256"]
    p0=S0["lab_momentum_edges_GeV"][-1]
    b0=p0/np.sqrt((E0-mp)*(E0+mp))
    bg=TabulatedRadiation(eos,S0["T_visible_GeV"])
    p_first=S0["lab_momentum_edges_GeV"][1]
    def gap(b):
        return window(p_first/b)[-2]-np.hypot(mp,p0/b)
    full_cell_cut=brentq(gap,b0,b0*1.02,xtol=2e-14)
    curves=[]
    for j in range(1,limit+1):
        endpoints=[read(work/f"q{q}"/f"step-{j}.json") for q in (32,64)]
        S,T=[m["snapshot"] for m in endpoints]
        assert all(m["snapshot_sha256"]==digest(m["snapshot"]) for m in endpoints)
        assert abs(S["T_visible_GeV"]/T["T_visible_GeV"]-1)<1e-14
        np.testing.assert_allclose(S["lab_momentum_edges_GeV"],
                                   T["lab_momentum_edges_GeV"][:len(S["lab_momentum_edges_GeV"])],
                                   rtol=1e-13,atol=0.)
        common=len(S["lab_momentum_edges_GeV"])
        crop=dict(T,lab_momentum_edges_GeV=T["lab_momentum_edges_GeV"][:common],
                  entropy_variable_scalar=np.asarray(T["entropy_variable_scalar"])[:,:common-1].tolist(),
                  entropy_variable_fermion=np.asarray(T["entropy_variable_fermion"])[:,:common-1].tolist())
        re=[]
        for label,k in (("q32",S),("q64",T),("q64 cropped to q32",crop)):
            Q=freeze(k,eos,10,16)
            H=Q["H_over_rest_decay_width"]*width
            P=np.asarray(k["lab_momentum_edges_GeV"])
            th=[np.asarray(k["entropy_variable_"+k0]) for k0 in ("scalar","fermion")]
            J=Collision(P)
            active=int(np.sum(J.J==0))
            rate=float(J.rate(th)[1,0,0]/J.c[1].vol[0]/H)
            if Q["widest_known_scalar_E_GeV"]<Q["local_scalar_energy_window_GeV"][0]:
                assert active==0 and rate==0.
            N=get_mom(P,th)[0][1,0]
            d=dict(label=label,scalar_edge_GeV=Q["widest_known_scalar_E_GeV"],
                           first_light_mean_occupation=float(N/J.c[1].vol[0]),
                           rest_scalar_window_GeV=Q["local_scalar_energy_window_GeV"],
                           native_cold_cell_event_count=active,
                           native_first_light_cell_number_rate_per_Hubble=rate,
                           independent_open_rest_gain_stimulated_inverse_per_Hubble=
                             Q["spontaneous_stimulated_and_inverse_at_rest_per_Hubble"])
            if label=="q64" and (j in (1,limit) or S["scale_factor"]/S0["scale_factor"]>b0):
                avg=independent_first_average(k,H,16)
                assert abs(rate/avg-1.)<2e-5
                d["independent_open_first_shell_averaged_prompt_rate_per_Hubble"]=avg
            re.append(d)
        if curves:
            dt=np.log(S["scale_factor"]/S0["scale_factor"]/curves[-1]["relative_a"])
            for x,y in zip(re,curves[-1]["controls"],strict=True):
                x["mean_occupation_change_per_log_a_over_preceding_interval"]=(x["first_light_mean_occupation"]-y["first_light_mean_occupation"])/dt
        curves.append(dict(step=j,T_GeV=S["T_visible_GeV"],relative_a=S["scale_factor"]/S0["scale_factor"],
                           checkpoints=[a["snapshot_sha256"] for a in endpoints],controls=re))
    return dict(scope="Matched selected coupled checkpoints, with independent frozen prompt rate "
                      "audit. Native represented operator requires both daughters in support; "
                      "the open frozen diagnostic uses a stationary parent surrogate.",
                q32_rest_recoil_cross_b=b0,q32_rest_recoil_cross_T_GeV=bg.state(b0)[0],
                q32_whole_first_shell_cut_b=full_cell_cut,
                q32_whole_first_shell_cut_T_GeV=bg.state(full_cell_cut)[0],frames=curves)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path);p.add_argument("work",type=Path)
    p.add_argument("--tapes",type=Path,default=Path("results"));p.add_argument("--limit",type=int,default=16)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.tapes,a.work,a.limit),indent=2,allow_nan=False))
