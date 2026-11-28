"""Measure known gross exterior creation without confusing it with exterior reuse."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np

from frw_exterior_restart_reference import import_restart


def read(path):
    raw=Path(path).read_bytes()
    return json.loads(gzip.decompress(raw) if str(path).endswith(".gz") else raw)


def ext(state):
    return import_restart(state["positive_exterior_restart"])


def run_case(parent,tape,E,threshold,old_energy,require_sequential_labels=True):
    A=read(parent)
    assert A["blocks_finished"]==A["blocks_total"]==32
    S=A["snapshot"]
    ext0=ext(S)
    scale0=S["scale_factor"]
    initial_upper=ext0.upper_bound(1.)
    rows=[]
    for i,r in enumerate(tape if isinstance(tape,list) else read(tape),start=1):
        if require_sequential_labels:
            assert r["blocks_finished"]==i
        x=r["snapshot"]; current=ext(x)
        a=x["scale_factor"]/scale0
        assert a>(rows[-1]["b_relative"] if rows else 1.)
        N0=ext0.array()
        C=current.array()
        k=N0.shape[1]
        # The original cohorts must persist with the same mass, side,
        # energy-at-injection and comoving labels across rebase.
        np.testing.assert_allclose(C[[0,5],:k]*(np.array([a,a**3]))[:,None],N0[[0,5]],
                                   rtol=3e-14,atol=0.)
        assert np.array_equal(C[[1,2,3,4,6],:k],N0[[1,2,3,4,6]])
        Nnew=C[:,k:]
        created_energy=Nnew[5]*Nnew[6]*a**3
        # Exact injection energy in original a^4 units is weighted by
        # injection expansion relative to original anchor, not current b.
        delta=created_energy @ (Nnew[0]*a)
        raw=current.injected*a**4-ext0.injected
        assert abs(raw-delta)<1e-12*max(raw,delta,ext0.injected)
        high=Nnew[3]==1
        y=sum(created_energy[high])
        if rows:
            assert y>=rows[-1]["cumulative_high_gross_created_GeV4"]-1e-13*max(y,1e-280)
        pool=initial_upper+y
        rows.append(dict(stage=i,b_relative=a,tracked_high_new_cohorts=int(high.sum()),
                         total_tracking_new_cohorts=Nnew.shape[1],
                         cumulative_high_gross_created_GeV4=float(y),
                         cumulative_high_gross_created_over_available=float(y/E),
                         all_initial_known_exterior_plus_created_high_energy_upper_GeV4=float(pool),
                         older_unrepresented_pressure_only_energy_upper_GeV4=old_energy,
                         one_use_pool_including_older_energy_GeV4=pool+old_energy,
                         ratio_of_that_one_use_pool_to_needed_gross_limit=float((pool+old_energy)/(E*threshold))))
    return rows


def run(pions,peripheral,exterior,p16,p32,t16,t32):
    P,G,X=read(pions),read(peripheral),read(exterior)
    rows=[]
    for a,b,c,s,t in zip(P["cases"],G["cases"],X["cases"],(p16,p32),(t16,t32),strict=True):
        q=a["q"]
        assert q==b["q"]==c["q"]
        E=a["controls"][-1]["new_available_energy_upper_GeV4"]
        R=b["max_gross_comoving_peripheral_input_fraction_for_test"]
        rows.append(dict(q=q,required_gross_peripheral_fraction_for_partial_test=R,
                         controls=run_case(s,t,E,R,c["old_exterior_rho_pressure_only_upper_GeV4"])))
    return dict(scope="Measures positive one-way exterior creation for the saved warm "
                      "continuation, and gives a one-use energy pool only under "
                      "separate no-reuse/no-reheating assumptions. Not gross exterior "
                      "reaction use, which is absent from this integrator. For older "
                      "unresolved tails only the separately conditioned energy upper is included.",
                cases=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("pions","peripheral","exterior","q16_initial","q32_initial","q16_tape","q32_tape"):
        p.add_argument(k)
    a=p.parse_args()
    print(json.dumps(run(*vars(a).values()),indent=2,allow_nan=False))
