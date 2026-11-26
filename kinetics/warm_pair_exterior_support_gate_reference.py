"""Sharper warm support gate when exterior scalars are also unoccupied."""
import json
from pathlib import Path

import numpy as np

from three_species_decay_onsager_reference import m,mp
from warm_scalar_pair_late_support_bound_reference import report as broad_bound


def scalar_support(P):
    assert P>0 and m>mp
    # At least one scalar's energy exceeds the retained edge if edge<m;
    # two retained scalars also cannot make a stable pair in this case.
    assert np.hypot(P,mp)<m
    B=m*m-mp*mp
    S=m*m+mp*mp
    pcrit=np.sqrt(B*B/(S+np.sqrt(S*S+3*B*B)))
    bcrit=max(1.,P/pcrit)
    assert abs(np.hypot(m,pcrit)-pcrit-np.hypot(mp,pcrit))<1e-15
    return pcrit,bcrit


def trials(P,bcrit):
    rng=np.random.default_rng(1534)
    waves=[]
    for b in (1.,bcrit*(1.+1e-8),2.,3.):
        p=P/b
        n=50000
        P1,P2=rng.random((2,n))*p
        c=rng.uniform(-1,1,n)
        E=np.hypot(m,P1)+np.hypot(m,P2)
        k=np.sqrt(P1*P1+P2*P2+2*P1*P2*c)
        s=E*E-k*k
        s=np.maximum(s,4*m*m)
        e_minus=(E-k*np.sqrt(1.-4*mp*mp/s))/2
        e_plus=E-e_minus
        edge=np.hypot(p,mp)
        assert np.all(e_plus>edge)
        assert np.all(e_minus>=np.hypot(m,p)-p)
        if b>bcrit:
            assert np.all(e_minus>edge)
        waves.append(dict(b=b,minimum_sampled_e_minus_over_edge=float(min(e_minus)/edge),
                          exact_lower_e_minus_over_edge=float((np.hypot(m,p)-p)/edge)))
    return waves


def report(base,late=None):
    P=base["support_upper_GeV"]
    pcrit,bcrit=scalar_support(P)
    ref=base["stable_fixed_count_Omega_h2_at_anchor"]
    tau=base["quadratic_loss_integral_without_Bose_enhancement"]
    need=base["required_normalized_quadratic_loss_integral"]
    active=1.-bcrit**-3
    result=dict(both_final_scalars_outside_for_b_greater_than=bcrit,
                equal_scalar_and_stable_support_pcrit_GeV=pcrit,
                T_at_anchor_GeV=base["temperature_GeV"],
                inverse_closed_for_entire_retained_incident_support=True,
                included_fraction_of_geometric_integral_where_one_stimulated_final_scalar_possible=active,
                conditional_cases=[dict(cap=F,loss_integral=tau*(1.+F*active),
                                        lower_Omega_stable=ref/(1.+tau*(1.+F*active)))
                                   for F in (0.,.001,1.,2.,10.,100.)],
                sampled_kinematics=trials(P,bcrit))
    if active>0.:
        result["necessary_retained_occupation_supremum_to_allow_target"]=(need/tau-1.)/active
    if late is not None:
        assert late["maximum_midpoint_inverse_over_forward"]==0.
        result["late_concrete_comparator"]=dict(steps=late["steps"],
             loss_fraction=1.-late["endpoint"]["inclusive_counts_relative_to_start"][1],
             endpoint_temperature_GeV=late["endpoint"]["photon_temperature_GeV"],
             no_inverse=True)
    return result


if __name__=="__main__":
    import argparse
    p=argparse.ArgumentParser()
    p.add_argument("general_bound",type=Path)
    p.add_argument("--late-result",type=Path)
    p.add_argument("--end-shape",type=Path)
    p.add_argument("--end-count-audit",type=Path)
    a=p.parse_args()
    base=json.loads(a.general_bound.read_text())
    q=None if a.late_result is None else json.loads(a.late_result.read_text())["rows"][-1]
    ans=report(base,q)
    assert (a.end_shape is None)==(a.end_count_audit is None)
    if a.end_shape:
        s=json.loads(a.end_shape.read_text())["snapshots"][-1]
        r=json.loads(a.end_count_audit.read_text())["stages"][-1]
        omega=r["retained_and_all_particle_tail_Omega_h2"][1]
        assert s["steps"]==r["n_steps"] and abs(s["T_visible_GeV"]/r["temperature_GeV"]-1.)<1e-12
        c=broad_bound(s,omega,r["all_tail_Omega_h2"][1])
        _,bcrit=scalar_support(c["support_upper_GeV"])
        assert bcrit==1.
        ans["terminal_pair_only_interval_from_end_anchor"]=dict(
            anchor_temperature_GeV=c["temperature_GeV"],
            loss_integral_upper=c["quadratic_loss_integral_without_Bose_enhancement"],
            conditional_stable_Omega_h2_lower_and_upper=[
                c["conditional_cases"][0]["final_stable_fixed_count_Omega_h2_lower"],omega],
            inverse_closed=True,all_forward_scalar_final_states_outside=True)
    ans["scope"]="Equal finite comoving scalar and stable support; no incident or stimulating particles outside the imposed support. Physical exterior occupancy or growing support invalidates the strengthening."
    print(json.dumps(ans,indent=2,allow_nan=False))
