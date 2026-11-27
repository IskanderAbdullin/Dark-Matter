"""Transfer the charged-pion windows into explicit ballistic distance bands."""
import argparse,json
from pathlib import Path
import numpy as np
from archived_collisionless_streaming_reference import Background
from kloe_pion_band_portal_reference import table,NB_PER_GEV2
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import P
from timelike_portal_onsager_coverage_reference import V,width,gprod,mpi
from warm_babar_dispersive_gap_reference import POINT,mu_kernel_lower
from warm_late_cold_pion_dispersive_recoil_bound_reference import phase_bound
from warm_late_cold_visible_offpole_dose_reference import log_tail
from warm_light_occupation_stage_audit_reference import read


def far_tail(bg,q,nt):
    from decay_scalar_ir_relaxation_reference import limits
    from three_species_decay_onsager_reference import M
    T=np.r_[np.geomspace(.05,.012,nt+1),.0001]
    A=np.array([bg.scale(t) for t in T]);p=q/A[:-1]
    u=np.log(T[1:])
    denom=(2*np.pi**2/45)*T[1:]**5*np.exp(bg.gs(u)+.5*bg.gr(u))*np.sqrt(np.pi**2/90)/MPL
    angle=p*p/(2*P*(P+np.sqrt((P-p)*(P+p))))
    B=(np.hypot(V*V,V*width)+V*V)/(2*(V*width)**2)
    coeff=2*12*(3*gprod*gprod/(8*np.pi))*B/(120*np.pi**4)
    hot=T[:-1]
    limit=np.maximum(4.,np.array([limits(pp)[1]+M for pp in p]))
    X=np.array([log_tail(a/b) for a,b in zip(limit,hot,strict=True)])
    r=np.exp(np.log(coeff)+6*np.log(hot)+X)
    return float(np.diff(np.log(A))@(r*angle*((-1/np.expm1(-mpi/hot))**2)/denom))


def run(eos,distance,gap,raw):
    bg=Background(eos)
    d=read(distance)
    G=read(gap)["bands"]
    gd,cs,err=table(raw)
    cs=(cs+2*err)/NB_PER_GEV2
    rows=[]
    for p in d["rows"][:3]:
        q=bg.T0*p["conservative_late_placement_cut_q_over_today_T"]
        H=p["archive_held_distance_selected_yield"]
        afds=[]
        for nt,nq,nx in ((32,30,24),(72,48,42)):
            obs=phase_bound(bg,gd-.005,gd+.005,nt,nq,nx,qcut_today=q)
            kloe=float(.01*cs@obs)
            extra=0.
            for a in G:
                s=np.linspace(*a["s_GeV2"],25)
                v=phase_bound(bg,s[:-1],s[1:],nt,nq,nx,qcut_today=q)
                mu=mu_kernel_lower(s[1:])/(4*np.pi**3)
                extra+=float(np.max(v/mu))*a["assumed_positive_dispersion_cap"]
            s=np.geomspace(.95,16,66)
            v=phase_bound(bg,s[:-1],s[1:],nt,nq,nx,qcut_today=q)
            other=float(np.diff(s)@(v*POINT/s[:-1]))+far_tail(bg,q,nt)
            afds.append(dict(T_panels=nt,KLOE_conditional_upper=kloe,low_gap_conditional_upper=extra,
                             high_unit_R_upper=other,pion_total_at_high_R_1_over_held=(kloe+extra+other)/H))
        L=10**p["controls"][-1]["log10_direct_lepton_yield_upper"]
        rows.append(dict(max_distance_Mpc=p["max_path_Mpc"],eligible_q_over_T0=q/bg.T0,
                         held_yield=H,lepton_upper=L,pion_choices=afds,
                         combined_conditional_upper_over_held_at_high_R_1=L/H+afds[-1]["pion_total_at_high_R_1_over_held"]))
    return dict(scope="Ballistic travel only; charged-pion vacuum spectrum with KLOE fixed bin plus twice marginal "
                      "inputs, separate conditional low-gap moments and a chosen upper spectral multiplier at high s. "
                      "No unmeasured source, scattering, structure or confidence-limit claim.",cases=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","distance","gap","kloe"):
        p.add_argument(k,type=Path)
    a=p.parse_args();print(json.dumps(run(a.eos,a.distance,a.gap,a.kloe),indent=2,allow_nan=False))
