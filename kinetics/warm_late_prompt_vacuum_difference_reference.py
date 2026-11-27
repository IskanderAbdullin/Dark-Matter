"""Compare activation against a mesh- and fixed-clock-matched vacuum pad."""

import argparse
import json
from pathlib import Path

import numpy as np

from frw_exterior_restart_reference import import_restart
from three_species_decay_onsager_reference import m,mp
from warm_history_lower_matched_audit_reference import measures
from warm_late_prompt_matched_vacuum_pad_reference import velocity_moment
from warm_light_occupation_stage_audit_reference import digest,read


def run(source,active,vacuum,limit):
    seed=read(source)
    S=seed["snapshot"]
    assert seed["snapshot_sha256"]==digest(S)
    old=import_restart(S["positive_exterior_restart"]).array()
    outside=np.array([old[5,old[1]==mass].sum() for mass in (mp,m)])
    shared=len(S["lab_momentum_edges_GeV"])-1
    rows=[]
    anchors=None
    for i in range(limit+1):
        A=read(active/f"step-{i}.json")
        B=read(vacuum/f"step-{i}.json")
        X,Y=A["snapshot"],B["snapshot"]
        assert A["hash"]==digest(X) and B["hash"]==digest(Y)
        assert A["time"]==B["time"]
        assert A["plan"]["eos_sha256"]==B["plan"]["eos_sha256"]
        assert A["plan"]["substeps"]==B["plan"]["substeps"]
        assert A["plan"]["artificial_empty_shell_filling"]==0.
        assert B["plan"]["source_state_sha256"]==seed["snapshot_sha256"]
        assert B["plan"]["matching_active_grid_sha256"]==A["plan"]["seed_sha256"]
        assert B["plan"]["exact_matching_clock_from_active_case"]
        assert X["T_visible_GeV"]==Y["T_visible_GeV"] and X["scale_factor"]==Y["scale_factor"]
        p=np.asarray(X["lab_momentum_edges_GeV"])
        np.testing.assert_array_equal(p,Y["lab_momentum_edges_GeV"])
        assert not X["positive_exterior_restart"]["cohorts"]
        assert Y["positive_exterior_restart"]==S["positive_exterior_restart"]
        U,V=measures(X),measures(Y)
        N=V[0].sum(axis=1)
        assert np.max(abs((U[0]-V[0]).sum(axis=1)-outside)/N)<3e-13
        vx=np.array(velocity_moment(p,[np.asarray(X["entropy_variable_"+k]) for k in ("scalar","fermion")]))
        vy=np.array(B["velocity2_and_p2_p4_by_species"])
        if anchors is None:
            anchors=(U.copy(),V.copy(),vx.copy(),vy.copy())
        U0,V0,Vx0,Vy0=anchors
        s=V0.sum(axis=-1)
        diff=U-V
        changed=(U-U0)-(V-V0)
        rate=A["native_signed_cut_flux"]
        tail=sum(A["excess_over_initial_energy"])
        rows.append(dict(time=A["time"],active_snapshot_sha256=A["hash"],vacuum_snapshot_sha256=B["hash"],
                         normalized_full_moment_L1=(abs(diff).sum(axis=-1)/s).tolist(),
                         normalized_original_support_L1=(abs(diff[:,:,:shared]).sum(axis=-1)/s).tolist(),
                         normalized_change_in_signed_original_support_L1=(abs(changed[:,:,:shared]).sum(axis=-1)/s).tolist(),
                         velocity2_p2_p4_relative_active_vs_vacuum=(vx/vy-1.).tolist(),
                         velocity_relative_increment_difference=(((vx-vy)-(Vx0-Vy0))/Vy0).tolist(),
                         active_and_vacuum_cut_excess=[tail,B["cut_excess_over_matched_clock_energy"]],
                         active_and_vacuum_signed_cut_rate=[rate,B["signed_cut_flux_over_matched_energy"]]))
    return dict(scope="Same finite mesh and fixed clock, same prompt-only operator. The old receipts "
                      "are restored only in the active run. These velocity moments are not a Boltzmann "
                      "hierarchy, structure-formation transfer function, or whole-history bound.",
                hard_GeV=A["plan"]["hard_GeV"],species=["scalar","light"],moment_order=["number","energy","pressure"],
                velocity_moment_order=["v_squared","p_squared","p_to_fourth"],frames=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("source","active","vacuum"):
        p.add_argument(k,type=Path)
    p.add_argument("--limit",type=int,default=6)
    a=p.parse_args()
    print(json.dumps(run(a.source,a.active,a.vacuum,a.limit),indent=2,allow_nan=False))
