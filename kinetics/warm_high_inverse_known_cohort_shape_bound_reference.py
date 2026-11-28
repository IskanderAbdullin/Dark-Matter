"""Shape-independent initial high inverse source bound with finite admitted support."""

import argparse
import json
from pathlib import Path

import numpy as np

from frw_exterior_restart_reference import import_restart
from late_radiation_background_reference import MPL,TabulatedRadiation
from three_species_decay_onsager_reference import M,m,mp
from factorized_prompt_return_graph_reference import C
from warm_cohort_two_packet_geometry_reference import highrate,mom,seed
from warm_light_occupation_stage_audit_reference import digest,read


def bound(s,eos,ratio):
    p=s["lab_momentum_edges_GeV"][-1]
    mm=np.array([mp,m]); deg=np.array([1.,4.])
    maximum=np.hypot(ratio*p,mm)
    threshold=np.hypot(M,1.)
    Emin=threshold-maximum[::-1]
    assert np.all(Emin>np.hypot(p,mm))  # both incoming particles must be from outside the retained q32 grid
    k=np.sqrt((Emin-mm)*(Emin+mm))
    ext=import_restart(s["positive_exterior_restart"])
    a,mass,_,side,edge,n,mean=ext.array()
    assert np.all(side==1)
    assert np.all(np.isin(mass,mm))
    count=[]
    for q,z in zip(mm,k,strict=True):
        mask=mass==q
        req=np.hypot(q,z/a[mask])
        lower=edge[mask]
        # E is integrated on the cohort's unknown positive distribution.
        # The rest-frame shape or correlations within it are not reconstructed.
        assert np.all(req>lower)
        share=np.clip((mean[mask]-lower)/(req-lower),0.,1.)
        count.append(float(n[mask]@share))
    integrals=np.asarray(count)*(2*np.pi**2/deg)/(k*Emin)
    rr=mom(s)[2].sum()
    H=np.sqrt((TabulatedRadiation(eos,s["T_visible_GeV"]).state(1.)[1]+rr)/3)/MPL
    # Omitting final heavy Pauli can only increase the rate. In the
    # independent fixed-E_H measure, dE_H dE_phi = dE_L dE_phi.
    # Drop the exact angular constraint, retaining the necessary individual
    # thresholds and upper bound on total E from finite admitted support.
    z=C*maximum.sum()*np.prod(integrals)/(rr*H)
    return dict(ratio_outer_to_q32_support=ratio,admitted_max_momentum_GeV=ratio*p,
                each_incoming_energy_necessary_GeV=Emin.tolist(),
                per_species_eligible_particle_number_upper_GeV3=count,
                normalized_initial_high_inverse_energy_upper=float(z))


def run(eos,D):
    R=read(D/"warm-light-occupation-next-q32-tape-2026-11-26.json.gz")[0]
    s=R["snapshot"]; assert R["snapshot_sha256"]==digest(s)
    out=[bound(s,eos,b) for b in (1.75,2.,2.15)]
    ref=out[1]["normalized_initial_high_inverse_energy_upper"]
    controls=[]
    for p in (.52,.58):
        x,_=seed(s,p)
        y=highrate(x,1.,eos)[-1][0]
        assert 0.<y<ref
        controls.append(dict(hard_packet_momentum_GeV=p,
                             initial_normalized_high_inverse_energy=y,upper_over_test=ref/y))
    return dict(scope="At the archived anchor only, specified bare inverse decay and known exterior "
                      "cohorts with arbitrary birth distributions of the same number and mean energy, "
                      "conditioned on admitted present support and collisionless transport before "
                      "activation. No untracked tails, further re-heating, recycling, or total heavy bound.",
                snapshots_sha256=R["snapshot_sha256"],controls=out,explicit_moment_matched_tests=controls)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("eos",type=Path)
    p.add_argument("--tapes",type=Path,default=Path("results"));x=p.parse_args()
    print(json.dumps(run(x.eos,x.tapes),indent=2,allow_nan=False))
