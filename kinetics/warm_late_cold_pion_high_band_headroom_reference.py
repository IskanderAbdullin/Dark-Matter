"""Translate an unspecified high-invariant charged-pion cap into cold-yield headroom."""
import argparse,json
from pathlib import Path
import numpy as np
from archived_collisionless_streaming_reference import Background
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import P
from timelike_portal_onsager_coverage_reference import V,width,gprod,mpi
from warm_babar_dispersive_gap_reference import POINT
from warm_late_cold_pion_dispersive_recoil_bound_reference import phase_bound
from warm_late_cold_visible_offpole_dose_reference import log_tail
from warm_light_occupation_stage_audit_reference import read


def tail(bg,nt):
    T=np.r_[np.geomspace(.05,.012,nt+1),.0001]
    A=np.array([bg.scale(t) for t in T])
    cut=.1*bg.T0/A[:-1]
    angle=cut**2/(2*P*(P+np.sqrt((P-cut)*(P+cut))))
    # A unit pion R(s) uses the same invariant weight as one massless
    # MB charged-lepton channel. Reuse its loose three-channel global
    # maximum, then multiply by the physical pion BE upper.
    bw=(np.hypot(V*V,V*width)+V*V)/(2*(V*width)**2)
    coeff=2*12*(3*gprod*gprod/(8*np.pi))*bw/(120*np.pi**4)
    raw=np.exp(np.log(coeff)+6*np.log(T[:-1])+np.array([log_tail(4/t) for t in T[:-1]]))
    BE=(-1/np.expm1(-mpi/T[:-1]))**2
    Te=T[1:];ln=np.log(Te)
    H=Te**2*np.sqrt(np.pi**2/90*np.exp(bg.gr(ln)))/MPL
    s=(2*np.pi**2/45)*Te**3*np.exp(bg.gs(ln))
    return float(np.diff(np.log(A))@(raw*BE*angle/(H*s)))


def run(eos,fold,box,kloe,gap,overlap):
    bg=Background(eos)
    y=read(fold)["held_equality_cold_yield"]
    u=read(box)["rows"][-1]["controls"][-1]["upper_fraction_of_held_yield"]
    k=read(kloe)["controls"][-1]
    kchoice=(k["BE"]+2*k["linear_published_marginal_variation"])/y
    lowergap=read(gap)["sum_of_conditional_outside_KLOE_cold_yield_upper"]/y
    K=read(overlap)["scenarios"][2]
    assert K["compatible"] and K["held_fraction_upper"]<lowergap
    results=[]
    for ns,nt,nq,nx in ((40,32,24,24),(100,72,42,42)):
        edges=np.geomspace(.95,16,ns+1)
        vals=phase_bound(bg,edges[:-1],edges[1:],nt,nq,nx)
        # R = sigma * s / POINT, so the largest sigma per unit cap
        # in this s interval occurs at its lower edge.
        Y=float(np.diff(edges)@(vals*POINT/edges[:-1]))
        Z=tail(bg,nt)
        results.append(dict(s_panels=ns,T_panels=nt,held_fraction_upper_per_unit_R=(Y+Z)/y,
                            tail_s_above_16_held_fraction_upper=Z/y,finite_part_yield_upper=Y))
    C=results[-1]["held_fraction_upper_per_unit_R"]
    baseline=u+kchoice+lowergap
    return dict(scope="Conditional high-band vacuum pion form-factor cap. KLOE fixed bins are a separate chosen "
                      "spectrum, not a continuous confidence upper. No full thermal hadronic bound.",
                high_s_min_GeV2=.95,lepton_box_upper=u,KLOE_fixed_bin_plus_two_marginal_choice=kchoice,
                uncovered_low_independent_conditional_upper=lowergap,
                alternative_joint_overlap_low_upper=K["held_fraction_upper"],rows=results,
                example_total_with_uniform_high_R_1=baseline+C,
                uniform_high_R_cap_sufficient_for_conditional_total_below_held=(1-baseline)/C)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","fold","box","kloe","gap","overlap"):
        p.add_argument(k,type=Path)
    a=p.parse_args();print(json.dumps(run(a.eos,a.fold,a.box,a.kloe,a.gap,a.overlap),indent=2,allow_nan=False))
