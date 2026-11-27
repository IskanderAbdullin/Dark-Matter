"""Verify published tail and invariants independently from the saved finite-parent states."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from joint_prompt_two_moment_evolution_reference import Collision
from three_species_decay_onsager_reference import M,P,width
from warm_corrected_late_prompt_cut_aligned_reference import excess
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom


def audit(raw, frozen):
    A=read(raw)
    S=A["snapshot"]
    assert digest(S)==A["hash"]
    base=read(frozen/"summary.json")
    assert base["initial_hash"]==A["hash"]
    p=np.asarray(S["lab_momentum_edges_GeV"])
    k=Collision(p)
    E=k.E
    P0=np.sqrt((E-M)*(E+M))
    # Direct quadrature measure: no cosmological rate normalization enters.
    measure=k.w*(2.*P*E*P0)/(width*M*M)
    def packet(s):
        return [np.asarray(s["entropy_variable_"+name]) for name in ("scalar","fermion")]
    def totals(s,h):
        n,_,e=get_mom(p,packet(s))
        return np.r_[n.sum(axis=1)+measure@h,e.sum()+(measure*E)@h]
    y0=np.zeros_like(E)
    if base["initial_parent"]=="equilibrium":
        v=packet(S)
        f=k.c[0].f(v[0],k.ph,k.I);g=k.c[1].f(v[1],k.ch,k.J)
        a=np.bincount(k.y,weights=k.pre*f*g,minlength=k.n)
        b=np.bincount(k.y,weights=k.pre*(1+f)*(1-g),minlength=k.n)
        y0=np.divide(a,a+b,out=np.zeros_like(a),where=a+b>0.)
    before=totals(S,y0)
    C=A["plan"]["normalized_cut_GeV"]
    rho=get_mom(p,packet(S))[2].sum()
    results=[]
    for j,row in enumerate(base["rows"],1):
        z=read(frozen/f"step-{j}.json")
        assert z["hash"]==digest(z["snapshot"]) and z["row"]==row
        y=np.asarray(z["parent_occupations"],dtype="<f8")
        assert len(y)==len(E) and np.all((y>=0.)&(y<=1.))
        assert np.array_equal(z["snapshot"]["lab_momentum_edges_GeV"],p)
        rest=totals(z["snapshot"],y)/before-1.
        assert max(abs(rest))<2e-12
        W=sum(excess(z["snapshot"],C))/rho
        assert abs(W/row["actual_W"]-1.)<2e-12
        assert row["parent_weighted_residual_over_binary64_rounding_scale"]<12.
        assert row["total_entropy_change_relative"]>=0.
        results.append(dict(time=row["time"],daughter_sha256=z["hash"],
                            parent_little_endian_float64_sha256=hashlib.sha256(y.tobytes()).hexdigest(),
                            conservation=rest.tolist(), independently_folded_W=W,
                            W_relative_shift=W/row["stationary_W"]-1.,
                            actual_and_stationary_flux=[row["actual_signed_flux"],row["stationary_signed_flux"]]))
    return dict(seed=A["hash"],initial_parent=base["initial_parent"],parent_native_energy_nodes=len(E),rows=results)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("raw",type=Path)
    p.add_argument("frozen",type=Path)
    a=p.parse_args()
    print(json.dumps(audit(a.raw,a.frozen),indent=2,allow_nan=False))
