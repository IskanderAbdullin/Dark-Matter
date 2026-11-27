"""Compare a matched vacuum pad with outgoing angles either closed or opened."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from joint_visible_heavy_adiabatic_reference import snapshot_functions
from three_species_decay_onsager_reference import m,mp
from warm_comoving_parent_chart_reference import ComovingCollision,common_rule
from warm_corrected_late_prompt_cut_aligned_reference import excess
from warm_finite_parent_comoving_expansion_reference import parent_equilibrium
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom,pad


def moments_in_pad(s,old,C):
    p=np.asarray(s["lab_momentum_edges_GeV"])[old:]
    nodes,weights=roots_legendre(32)
    dp=np.diff(p)[:,None]
    q=p[:-1,None]+dp*(nodes+1.)/2.
    B,_,F,_=snapshot_functions(s)
    val=[]
    for mass,deg,fun in ((mp,1,B),(m,4,F)):
        en=np.hypot(q,mass)
        assert np.all(en>C)
        n=deg/(4*np.pi**2)*dp*q**2*weights*fun(en)
        val.append([n.sum(),np.sum(n*en),np.sum(n*(en-C))])
    return np.array(val)


def run(anchor,base,opened):
    r=read(anchor)
    S=r["snapshot"]
    assert digest(S)==r["hash"]
    p0=np.asarray(S["lab_momentum_edges_GeV"])
    z=[read(t/"step-1.json") for t in (base,opened)]
    b,c=[t["row"] for t in z]
    assert not b["outer_prompt_is_open"] and c["outer_prompt_is_open"]
    assert b["time_J"]==c["time_J"] and b["step_H"]==c["step_H"] and b["c"]==c["c"]
    assert b["qsha256"]==c["qsha256"] and b["wsha256"]==c["wsha256"]
    C=r["plan"]["normalized_cut_GeV"]
    rho=get_mom(p0,pad(S,p0))[2].sum()
    old=len(p0)-1
    vv=[]
    for file,a in zip((base,opened),z,strict=True):
        P=a["daughter_snapshot"]
        row=a["row"]
        assert digest(P)==row["daughter_hash"] and row["initial_daughter_hash"]==r["hash"]
        assert row==read(file/"summary.json")["rows"][-1]
        assert np.allclose(np.asarray(P["lab_momentum_edges_GeV"])[:len(p0)]*row["c"],p0,rtol=3e-15)
        aa=moments_in_pad(P,old,C)
        np.testing.assert_allclose(aa[:,0]*row["c"]**3,row["newly_resolved_number_by_species_a3"],rtol=4e-11,atol=0.)
        np.testing.assert_allclose(aa[:,1:]/rho,np.array([row["newly_resolved_energy_by_species_over_rho"],
                            row["newly_resolved_excess_by_species_over_rho"]]).T,rtol=4e-11,atol=0.)
        x=sum(excess(P,C))/rho
        assert abs(x/row["W_over_initial_daughter_rho"]-1.)<2e-14
        vv.append(dict(total_W_over_started_rho=x,pad_N_E_W_by_species=aa.tolist(),daughter_hash=row["daughter_hash"],
                       energy_defect=row["relative_energy_defect"],root_residual=row["solver"]["residual"]))
    assert np.all(np.asarray(vv[0]["pad_N_E_W_by_species"])==0.) and np.all(np.asarray(vv[1]["pad_N_E_W_by_species"])>0.)
    cap=float(np.asarray(z[1]["daughter_snapshot"]["lab_momentum_edges_GeV"])[-1]*c["c"])
    q=np.r_[p0,np.geomspace(p0[-1],cap,4)[1:]]
    # Recover the initial occupied phase independently and record its identity.
    # This pad test uses the union at initial, midpoint and final scales.
    scales=np.exp(np.linspace(0.,c["step_H"],3))
    Q,w=common_rule(q,scales,read(opened/"summary.json")["radial_order"])
    def sha(v):
        return hashlib.sha256(np.asarray(v,dtype="<f8").tobytes()).hexdigest()
    if sha(Q)!=c["qsha256"]:
        # The endpoint redshift can round an exact decimal cap by one ulp.
        # Try the rounded cap used by the caller without refitting any phase.
        cap=round(cap,12)
        q=np.r_[p0,np.geomspace(p0[-1],cap,4)[1:]]
        Q,w=common_rule(q,scales,read(opened/"summary.json")["radial_order"])
    assert sha(Q)==c["qsha256"] and sha(w)==c["wsha256"]
    initial=pad(S,q)
    assert np.all(get_mom(q,initial)[0][:,old:]==0.)
    eq=parent_equilibrium(ComovingCollision(q,1.,Q,w,retained_shells=old),initial)
    summarize=lambda key:np.asarray(vv[1]["pad_N_E_W_by_species"])[:,key].sum()/rho
    return dict(hard_GeV=r["plan"]["hard_GeV"],phase_sha256=r["hash"],physical_cut_GeV=C,
                padded_comoving_edge_GeV=cap,initial_parent_sha256=sha(eq),
                opened_and_closed_charts=[vv[1],vv[0]],opened_pad_full_energy_over_rho=summarize(1),
                opened_pad_excess_over_rho=summarize(2),
                total_W_relative_change=vv[1]["total_W_over_started_rho"]/vv[0]["total_W_over_started_rho"]-1.)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("anchor","closed","opened"):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.anchor,a.closed,a.opened),indent=2,allow_nan=False))
