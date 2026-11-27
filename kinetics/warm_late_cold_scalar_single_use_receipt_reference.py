"""Conditional one-use cold-scalar supply from recorded exterior cohorts only."""

import argparse
import json
from pathlib import Path

import numpy as np

from frw_exterior_restart_reference import import_restart
from joint_prompt_two_moment_evolution_reference import kin
from three_species_decay_onsager_reference import Es,M,P,m,mp
from warm_heavy_eligible_exterior_receipt_reference import tail
from warm_late_cold_scalar_prompt_partition_reference import counts,recoil_gate
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom


def angular_vacuum_bound(q):
    ec=np.hypot(mp,q)
    y0=np.arcsinh(P/mp)
    yc=np.arcsinh(q/mp)
    assert 0.<yc<y0
    assert Es>ec
    ystar=np.arccosh(Es/ec)

    def share(y):
        return np.maximum(ec-mp*np.cosh(y-y0),0.)/(2*P*np.sinh(y))

    assert y0-yc<ystar<y0+yc
    optimum=float(share(ystar))
    grid=np.linspace(y0-yc,y0+yc,4001)
    rates=share(grid)
    assert optimum>0. and np.max(rates)<=optimum*(1+2e-13)
    assert abs(share(y0-yc))<3e-15*optimum and abs(share(y0+yc))<3e-15*optimum
    exact_parent_gate=float(recoil_gate(np.array([q]),1.)[0][0])
    assert abs(M*np.cosh(y0-yc)/exact_parent_gate-1.)<6e-14
    return optimum, M*np.cosh(ystar)


def eligible(E,e):
    birth_ratio,mass,_,side,bound,N,mean=e
    assert np.all(E>mass)
    E_birth=np.hypot(mass,np.sqrt((E-mass)*(E+mass))/birth_ratio)
    counts=[]
    for s,b,n,mu,x in zip(side,bound,N,mean,E_birth,strict=True):
        if s==0 and b>=x:
            # Safe bound if a finite lower-side cohort can cross threshold.
            counts.append(n)
        else:
            counts.append(tail(s,b,n,mu,x)[0])
    v=np.asarray(counts)
    assert np.all((v>=0)&(v<=N))
    return [float(v[mass==mp].sum()),float(v[mass==m].sum())]


def run(paths):
    out=[]
    for raw in paths:
        x=read(raw)
        S=x["snapshot"]
        assert x["snapshot_sha256"]==digest(S)
        p=np.asarray(S["lab_momentum_edges_GeV"])
        e=import_restart(S["positive_exterior_restart"]).array()
        scalar_ab=np.asarray(S['entropy_variable_scalar'])
        N0=get_mom(p,[S["entropy_variable_"+k] for k in ("scalar","fermion")])[0][0].sum()
        r=[]
        for c in (.03,.1,.3):
            q=c*S["T_visible_GeV"]
            cold=float(counts(S,np.array([q]),24)[0])
            assert 0.<cold<N0
            v,Hstar=angular_vacuum_bound(q)
            gate,_,open_cap,_=recoil_gate(np.array([q]),p[-1])
            assert gate[0]>open_cap
            # Every heavy parent in this band consumes at least one
            # eligible exterior input; if only one is exterior, its
            # other partner is bounded by the represented support.
            E_ret_other=np.where(e[1]==mp,np.hypot(m,p[-1]),np.hypot(mp,p[-1]))
            all_reentry=eligible(gate[0]/2,e)
            only_one=eligible(gate[0]-E_ret_other,e)
            # Within a reconstructed entropy-affine cell, f is
            # monotone, so both one-sided endpoints suffice.
            included=p[:-1]<q
            extent=kin(np.minimum(p[1:][included],q),mp)-kin(p[:-1][included],mp)
            step=np.diff(kin(p,mp))[included]
            assert np.all((extent>0.)&(extent<=step))
            a,b=scalar_ab[:,included]
            z=float(np.min(np.minimum(a,a+b*extent/step)))
            ceiling=float(np.exp(-z)/(-np.expm1(-z)))
            cap_B=min(1.,2.*(1.+ceiling)*v)
            r.append(dict(q_over_T=c,E_parent_gate_GeV=float(gate[0]),
                parent_at_largest_vacuum_probability_GeV=Hstar,max_unpolarized_vacuum_angle_share=v,
                archived_cold_BE_max_f=ceiling,conditional_max_probability_if_f_light_le_half=cap_B,
                eligible_tracked_one_use_species_counts_GeV3=all_reentry,
                eligible_if_second_input_represented_GeV3=only_one,
                geometric_one_decay_upper_over_initial_cold_count=sum(all_reentry)/cold,
                vacuum_decay_upper_over_initial_cold_count=v*sum(all_reentry)/cold,
                half_FD_frozen_cold_Bose_upper_over_initial_cold_count=cap_B*sum(all_reentry)/cold,
                represented_scalar_N_GeV3=float(N0),represented_held_cold_N_GeV3=cold))
        out.append(dict(hash=x["snapshot_sha256"],tracked_cohorts=e.shape[1],controls=r))
    return dict(scope="Single use of stored positive exterior inputs at the frozen archive; no untracked source, "
                      "repeated circulation, acceleration or transport. Vacuum decay and the optional "
                      "global light occupation cap refer to open, unpolarized prompt decays.",cases=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("archives",nargs="+",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.archives),indent=2,allow_nan=False))
