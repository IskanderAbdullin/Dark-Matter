"""Ballistic-only diagnostics for the late soft scalar band and held rest mixture."""

import argparse
import json
from pathlib import Path

import numpy as np

from archived_collisionless_streaming_reference import Background
from archived_scalar_pair_reaction_reference import points
from frw_exterior_restart_reference import import_restart
from joint_visible_heavy_adiabatic_reference import snapshot_functions
from three_species_decay_onsager_reference import m,mp
from warm_late_cold_scalar_prompt_partition_reference import counts
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom


def phase_split(bg, edges, ai, mass, kmax, target):
    # Split based on available phase advance rather than re-fitting occupations.
    # A tiny substitute at p=0 is used only to estimate the split, not as
    # an integration floor.
    d=bg.streams(np.maximum(edges,1e-90)*ai,mass,ai,bg.aeq,4)
    parts=np.maximum(1,np.ceil(np.diff(d)*kmax/target).astype(int))
    return np.r_[*(np.linspace(a,b,int(n),endpoint=False)
                   for a,b,n in zip(edges[:-1],edges[1:],parts,strict=True)),edges[-1]]


def run(table,files):
    bg=Background(table)
    ks=np.array([.003,.03,.3,1.,3.,10.,30.,100.])
    arr=[]
    for path in files:
        d=read(path)
        S=d["snapshot"]
        assert digest(S)==d["snapshot_sha256"]
        T=S["T_visible_GeV"]
        ai=float(bg.scale(T))
        eq=bg.aeq
        # Equality in this fixed background follows the post-annihilation
        # constant visible entropy branch. Restricted below last EOS knot.
        Teq=bg.T0/eq
        assert np.log(Teq)<bg.u[0]
        p=np.asarray(S["lab_momentum_edges_GeV"])
        q_archive=.1*T
        q_equal_archive=.1*bg.T0/ai
        assert q_archive<q_equal_archive<p[-1]
        N,_,_=get_mom(p,[S["entropy_variable_"+name] for name in ("scalar","fermion")])
        N=N.sum(axis=1)
        e=import_restart(S["positive_exterior_restart"]).array()
        n_ext=np.asarray([e[5,e[1]==d].sum() for d in (mp,m)])
        masses=np.array([mp,m])
        fraction=lambda n: float(mp*n[0]/(masses@n))
        fs=fraction(N)
        # The unknown momenta do not enter the known exterior rest count.
        fs_total=fraction(N+n_ext)
        f,_,g,_=snapshot_functions(S)
        controls=[]
        for r,order,phase in ((12,4,6.),(20,7,3.)):
            obs=[]
            for idx,(d,gamma,F) in enumerate(((mp,1,f),(m,4,g))):
                edges=np.unique(np.r_[p,q_archive,q_equal_archive]) if idx==0 else p
                edges=phase_split(bg,edges,ai,d,ks[-1],phase)
                pp,_,ww=points(edges,d,F,r)
                ww*=gamma
                assert abs(ww.sum()/N[idx]-1.)<4e-11
                distance=bg.streams(pp*ai,d,ai,eq,order)
                sinc=np.sinc(ks[:,None]*distance/np.pi)
                B=(sinc@ww)/ww.sum()
                mean=float(ww@distance/ww.sum())
                small=[]
                if idx==0:
                    for label,th in (("at_anchor",q_archive),("at_equality",q_equal_archive)):
                        col=pp<th
                        frac=ww[col].sum()/ww.sum()
                        n_di=counts(S,np.array([th]),24)[0]
                        assert abs(frac*ww.sum()/n_di-1.)<4e-12
                        BB=sinc[:,col]@ww[col]/ww[col].sum()
                        small.append(dict(cold_at=label,scalar_count_fraction=float(frac),
                              mean_comoving_distance_Mpc=float(ww[col]@distance[col]/ww[col].sum()),
                              contribution_to_scalar_distance_first_moment=float(ww[col]@distance[col]/(ww@distance)),
                              scalar_cold_only_ballistic_monopole=BB.tolist()))
                obs.append(dict(mean_path_Mpc=mean,ballistic_monopole=B.tolist(),phase_split_shells=len(edges)-1,bands=small))
            cast=np.asarray(obs[1]["ballistic_monopole"])
            mixed=cast+fs*(np.asarray(obs[0]["ballistic_monopole"])-cast)
            delta=mixed-cast
            assert max(abs(delta))<=2*fs
            controls.append(dict(momentum_order=r,expansion_order=order,max_target_phase_per_panel=phase,species=obs,
                         full_minus_stable_ballistic_monopole=delta.tolist(),
                         largest_magnitude_absolute_delta=float(max(abs(delta)))))
        np.testing.assert_allclose([a["mean_path_Mpc"] for a in controls[0]["species"]],
                                   [a["mean_path_Mpc"] for a in controls[1]["species"]],rtol=4e-8)
        for j in range(2):
            np.testing.assert_allclose(controls[0]["species"][j]["ballistic_monopole"],
               controls[1]["species"][j]["ballistic_monopole"],rtol=0,atol=4e-9)
        # Cold oscillatory transforms need separate refinement checks.
        for j in range(2):
            np.testing.assert_allclose(controls[0]["species"][0]["bands"][j]["scalar_cold_only_ballistic_monopole"],
              controls[1]["species"][0]["bands"][j]["scalar_cold_only_ballistic_monopole"],rtol=0,atol=2e-9)
        arr.append(dict(archive=digest(S),T_visible_GeV=T,scale_factor=ai,
               momentum_cut_for_equal_p_over_T_divided_by_anchor_cut=float(q_equal_archive/q_archive),
               scalar_rest_fraction_of_represented_two_species=fs,
               scalar_rest_fraction_with_recorded_exterior_counts=fs_total,
               relative_change_with_known_exterior=fs_total/fs-1.,rows=controls))
    return dict(scope="Frozen number-preserving collisionless propagation only; for initial perturbations common "
               "to both species, spherical monopole uses j0(k*travelled_distance), without gravity. "
               "Not a transfer spectrum, exclusion, late production estimate, or limit on untracked populations.",
               equality_scale=bg.aeq,k_per_Mpc=ks.tolist(),cases=arr)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("archives",nargs="+",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.archives),indent=2,allow_nan=False))
