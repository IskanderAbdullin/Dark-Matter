"""Evaluate tracked and newly outgoing excess energy on the common pair chart."""

import argparse
import json
from pathlib import Path

import numpy as np

from frw_exterior_restart_reference import import_restart
from joint_visible_heavy_adiabatic_reference import snapshot_functions
from late_radiation_background_reference import MPL,TabulatedRadiation
from three_species_decay_onsager_reference import M
from warm_light_occupation_stage_audit_reference import digest,read
from warm_pair_common_fermion_chart_reference import rate
from warm_prompt_reentry_shape_audit_reference import get_mom


def run(eos,pressure,archive):
    q=read(archive)
    S=q["snapshot"]
    assert q["snapshot_sha256"]==digest(S)
    p=np.asarray(S["lab_momentum_edges_GeV"])
    B,bs,F,fs=snapshot_functions(S)
    hi=np.hypot(M,read(pressure)["heavy_momentum_cut_at_late_anchor_GeV"])
    cuts=(.10,.30,.50,hi)
    th=[S["entropy_variable_"+k] for k in ("scalar","fermion")]
    rho=get_mom(p,th)[2].sum()
    H=np.sqrt((TabulatedRadiation(eos,S["T_visible_GeV"]).state(1.)[1]+rho+
               import_restart(S["positive_exterior_restart"]).energy(1.))/3)/MPL
    rows=[]
    for order in ((3,7,6,2),(4,10,10,4),(6,14,14,6)):
        out={}
        _,_,gross,defect,lim=rate(p,F,B,order,excess_cut=cuts,excess_result=out)
        v=np.array(out.pop("rate_GeV5"))/(rho*H)
        beyond=np.array(cuts)>=max(bs[-1],fs[-1])
        np.testing.assert_allclose(v[beyond,2],v[beyond,3],rtol=4e-13,atol=1e-200)
        rows.append(dict(quadrature=order,**out,rate_over_retained_rhoH=v.tolist(),
                         gross_pair_energy_rates_over_retained_rhoH=(gross[:,1]/(rho*H)).tolist(),
                         conservative_pair_ledger_defect=defect,max_empty_fermion_p_GeV=lim))
    return dict(scope="Frozen resolved input phases, empty incoming exterior with permitted outgoing states. "
                      "Exact per-event convex energy moment, not a refitted archive moment or time integral.",
                snapshot_sha256=q["snapshot_sha256"],controls=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","pressure","archive"):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.pressure,a.archive),indent=2,allow_nan=False))
