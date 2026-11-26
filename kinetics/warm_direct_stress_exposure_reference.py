"""Integrate bare high-band inverse exposure in the isolated fixed-clock stress."""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid

from joint_prompt_two_moment_evolution_reference import Collision
from joint_resonant_frw_split_reference import pack
from late_archived_shape_continuation_reference import collide
from late_radiation_background_reference import MPL, TabulatedRadiation
from warm_cohort_two_packet_geometry_reference import highrate, mom, seed
from warm_direct_prompt_stress_reference import DirectPrompt
from warm_light_occupation_stage_audit_reference import digest, read


def go(eos, record, hard, sub, terminal=1e-5):
    S = record["snapshot"]
    assert digest(S) == record["snapshot_sha256"]
    A, ledger = seed(S, hard)
    p = np.asarray(A["lab_momentum_edges_GeV"])
    keys = ("entropy_variable_scalar","entropy_variable_fermion")
    U = pack([A[k] for k in keys])
    n, _, E = mom(A)
    T = S["T_visible_GeV"]
    H = np.sqrt((TabulatedRadiation(eos,T).state(1.)[1] + E.sum())/3)/MPL
    base = n.sum(axis=1)
    model = DirectPrompt(Collision(p),H,[],True)
    # Early almost constant decade is included explicitly. Each log
    # decade has sub steps, with the last time equal across the ladder.
    grid = np.r_[0., np.geomspace(1e-10,terminal,5*sub+1)]
    rates = [highrate(A,1.,eos)[-1]]
    ledgers = []
    for j, (t0, t1) in enumerate(zip(grid[:-1],grid[1:],strict=True),1):
        U, res, nacc, nrej = collide(model,U,t1-t0,T,continuation=True,balanced=True)
        st = dict(A, **{k:f.tolist() for k,f in zip(keys,model.arrays(U),strict=True)})
        N, _, Q = mom(st)
        defect = [float(np.max(abs(N.sum(axis=1)/base-1.))),float(abs(Q.sum()/E.sum()-1.))]
        assert max(defect) < 4e-12
        rates.append(highrate(st,1.,eos)[-1])
        ledgers.append(dict(residual=res, accepted=nacc, rejected=nrej, moment_error=defect))
        if j == 1 or (j-1) % sub == 0:
            print(json.dumps(dict(hard=hard, subdivision=sub, time=float(t1), stage=j)),
                  file=sys.stderr, flush=True)
    exposure = cumulative_trapezoid(np.asarray(rates),grid,axis=0,initial=0.)
    half = np.maximum(abs(exposure[-1]),1e-280)
    return dict(scope="Isolated direct-occupation fixed-clock prompt. Cumulative quadrature estimate, not a bound.",
                hard=hard, samples_per_decade=sub, source_sha256=digest(S), seed=ledger,
                H_times=grid.tolist(), high_rate_rhoH=rates, integrated_energy_fraction=exposure.tolist(),
                exposure_by_1e_minus_7_over_final=(exposure[np.argmin(abs(grid-1e-7))]/half).tolist(),
                residual_ledgers=ledgers)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("--tapes",type=Path,default=Path("results"))
    p.add_argument("--hard",type=float,required=True)
    p.add_argument("--steps-per-decade",type=int,required=True)
    a=p.parse_args()
    s=read(a.tapes/"warm-light-occupation-next-q32-tape-2026-11-26.json.gz")[0]
    print(json.dumps(go(a.eos,s,a.hard,a.steps_per_decade),indent=2,allow_nan=False))
