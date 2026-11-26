"""Check a first-order extrapolation as an admissible moment diagnostic."""

import argparse
import json
from pathlib import Path

import numpy as np

from frw_exterior_restart_reference import import_restart
from joint_prompt_two_moment_evolution_reference import Cells
from three_species_decay_onsager_reference import m,mp
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom


def run(D):
    src=read(D/"warm-light-occupation-next-q32-tape-2026-11-26.json.gz")[0]["snapshot"]
    T0=src["scale_factor"]
    P=np.asarray(src["lab_momentum_edges_GeV"])
    first=get_mom(P,[src["entropy_variable_scalar"],src["entropy_variable_fermion"]])
    n0=first[0].sum(axis=1)*T0**3
    rho0=first[2].sum()*T0**3
    x0=len(src["positive_exterior_restart"]["cohorts"])
    state=[]; additions=[]
    for n in (2,4,8,16):
        a=read(D/f"warm-matched-explicit-support-q32-{n}-blocks-2026-11-26.json.gz")
        S=a["snapshot"]
        assert a["snapshot_sha256"]==digest(S) and a["blocks_total"]==a["blocks_finished"]==n
        p=np.asarray(S["lab_momentum_edges_GeV"]); A=S["scale_factor"]
        assert abs(A/T0/1.001-1)<4e-14
        np.testing.assert_allclose(p*1.001,P,rtol=5e-14,atol=0.)
        Q=np.asarray(get_mom(p,[S["entropy_variable_scalar"],S["entropy_variable_fermion"]])[:2])*A**3
        ext=import_restart(S["positive_exterior_restart"]).array()[:,x0:]
        state.append(Q)
        additions.append(np.array([ext[5,ext[1]==c].sum() for c in (mp,m)])*A**3)
    correction=[2*b-a for a,b in zip(state[:-1],state[1:],strict=True)]
    outside=[2*b-a for a,b in zip(additions[:-1],additions[1:],strict=True)]
    final=T0*1.001
    cases=[]
    for moment,ext in zip(correction,outside,strict=True):
        assert np.all(moment>0.) and np.all(ext>0.)
        assert abs((moment[0].sum()+ext.sum())/n0.sum()-1)<2e-12
        err=[]
        for j,d in enumerate((mp,m)):
            mesh=Cells(P/1.001,d,(1,4)[j])
            coeff,e,_=mesh.reconstruct(moment[:,j]/final**3)
            err.append(e)
        cases.append(dict(shell_reconstruction_error=err,
                          combined_number_residual=float(abs((moment[0].sum()+ext.sum())/n0.sum()-1))))
    def diff(a,b):
        c=a-b
        return [*(np.sum(abs(c[0]),axis=1)/n0),*(np.sum(abs(c[1]),axis=1)/rho0)]
    exact=[diff(a,b) for a,b in zip(state[:-1],state[1:],strict=True)]
    rich=[diff(a,b) for a,b in zip(correction[:-1],correction[1:],strict=True)]
    assert np.max(np.asarray(rich[-1])/exact[-1])<.2
    return dict(scope="Same finite q32 mesh and short clock. Richardson combinations are diagnostics, "
                      "not accepted steps or a rigorous error bound. Receipt counts are algebraic "
                      "totals only, not extrapolated exterior phase distributions.",
                metric_order=["scalar_number","light_number","scalar_kinetic_energy","light_kinetic_energy"],
                raw_successive_L1=exact,extrapolated_successive_L1=rich,
                latest_defect_reduction=(np.asarray(exact[-1])/rich[-1]).tolist(),
                tests=cases)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tapes",type=Path,default=Path("results"))
    print(json.dumps(run(p.parse_args().tapes),indent=2,allow_nan=False))
