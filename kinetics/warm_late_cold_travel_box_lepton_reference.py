"""One-sided direct-lepton checks at .5, 1 and 2 Mpc reference travel."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from archived_collisionless_streaming_reference import Background
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import mp
from warm_late_cold_distance_source_reference import analytic_late_q_enclosure,q_limit
from warm_late_cold_visible_recoil_box_upper_reference import log_slope_minimum,rate
from warm_light_occupation_stage_audit_reference import read


def sufficient_q(bg, earliest, D):
    ak=np.exp(bg.v[0]);ae=bg.aeq
    u=np.r_[np.log(.05), bg.u[(bg.u<np.log(.05)) & (bg.u>bg.u[0])][::-1],bg.u[0]]
    assert np.all(np.diff(u)<0.)
    A=np.array([bg.scale(t) for t in np.exp(u)])
    assert abs(A[0]/earliest-1)<2e-14 and abs(A[-1]/ak-1)<2e-14
    gs,gr=bg.gs(u),bg.gr(u)
    # PCHIP stays between its tabulated endpoints on each knot cell.
    early_rhs=np.sqrt(bg.ogamma*.5*np.exp(np.minimum(gr[:-1],gr[1:])+
                 4/3*(np.log(bg.gs0)-np.maximum(gs[:-1],gs[1:])))+bg.om*A[:-1]+bg.ol*A[:-1]**4)
    a=np.geomspace(ak,ae,2049)
    rhs=np.sqrt(bg.or0+bg.om*a[:-1]+bg.ol*a[:-1]**4)
    def path_upper(x):
        k=x*bg.T0/mp
        return bg.D0*k*(np.diff(np.arcsinh(A/k))@(1/early_rhs)+
                         np.diff(np.arcsinh(a/k))@(1/rhs))
    R=brentq(lambda x:path_upper(x)/D-1,1e-8,10.,xtol=1e-14)*(1-2e-7)
    assert path_upper(R)<D
    return R*bg.T0


def run(eos,exact):
    bg=Background(eos)
    assert bg.ol>0. and bg.om>0.
    slopes=[log_slope_minimum(bg.gs,3),log_slope_minimum(bg.gr,4)]
    O=read(exact)
    rows=[]
    for a in O["rows"][1:]:
        D=a["max_actual_path_Mpc"]
        top,low=.05,.012
        ve0=bg.scale(top);velate=bg.scale(.0001)
        lower=sufficient_q(bg,ve0,D)
        upper=analytic_late_q_enclosure(bg,velate,D,1024)
        q0,_=q_limit(bg,ve0,D)
        q1,_=q_limit(bg,velate,D)
        assert lower<q0<q1<upper
        Y=a["held_distance_selected_yield"]
        controls=[]
        for ntda,nq,nr,nx in ((48,36,42,28),(128,64,76,48)):
            T=np.geomspace(top,low,ntda+1)
            A=np.array([bg.scale(t) for t in T])
            s=2*np.pi**2/45*T**3*np.exp(bg.gs(np.log(T)))
            H=np.sqrt(np.pi**2/90*np.exp(bg.gr(np.log(T))))*T**2/MPL
            q_hi=np.array([sum(rate(t,upper/b,nq,nr,nx)[:2]) for t,b in zip(T[:-1],A[:-1],strict=True)])
            q_lo=np.array([rate(t,lower/b,nq,nr,nx)[2] for t,b in zip(T[1:],A[1:],strict=True)])
            dt=np.diff(np.log(A))
            ylo=float(dt@(q_lo/(s[:-1]*H[:-1])))
            yhi=float(dt@(q_hi/(s[1:]*H[1:])))
            expected=a["controls"][-1]["lepton_FD"]
            assert ylo<expected<yhi
            controls.append(dict(time_panels=ntda,Q_base=nq,R_panels=nr,angle_panels=nx,
               FD_lower_over_held=ylo/Y,MB_upper_over_held=yhi/Y))
        rows.append(dict(max_path_Mpc=D,known_sufficient_q_over_T0=float(lower/bg.T0),
                         known_necessary_q_over_T0=float(upper/bg.T0),
                         numerical_lepton_FD_over_held=a["controls"][-1]["lepton_FD"]/Y,controls=controls))
    return dict(monotone_EOS_slopes=slopes,scope="12–50 MeV direct-only lepton and unpolarized vacuum decay. Independent geodesic lower "
                      "and upper integrals bound the birth cut and kinematic boxes enclose the yield. "
                      "Not a complete hadronic or perturbation result.",cases=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path);p.add_argument("exact",type=Path)
    a=p.parse_args();print(json.dumps(run(a.eos,a.exact),indent=2,allow_nan=False))
