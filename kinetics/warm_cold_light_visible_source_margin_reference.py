"""Selected visible positive-rate margins for the frozen near-rest light control."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import expit,kve

from timelike_portal_onsager_coverage_reference import mpi
from three_species_decay_onsager_reference import M,m,width
from warm_late_visible_recycling_envelope_reference import COUP2,LEPTONS,VECTOR_MASS,VECTOR_GAMMA
from warm_light_occupation_stage_audit_reference import read


def bath(d,T):
    """Maxwell number and mean energy, per charge (two spin states)."""
    x=d/T
    return d*d*T*np.exp(-x)*kve(2,x)/np.pi**2, 3*T+d*kve(1,x)/kve(2,x)


def upper(d,T,Emin,Emax,Kpi):
    """Uniform momentum-mode injection bound on the specified energy band."""
    n,mean=bath(d,T)
    hot,_=bath(d,2*T)
    u=VECTOR_MASS**2
    w=VECTOR_MASS*VECTOR_GAMMA
    s=np.hypot(u,w)  # maximum of s/((s-u)^2+w^2) over s >= 0
    high=3*COUP2*s/(4*np.pi*((s-u)**2+w*w))
    low=3*COUP2/(4*np.pi*(u-1)**2)
    # Angle average of s = 2d^2+2(E E2-p p2 cos(theta)), after extending
    # the positive subpole integrand to all partner directions.
    up_low=low*n*np.exp(-Emin/T)*(2*d*d+2*Emax*mean)
    up_high=high*hot*np.exp(-(1+Emin)/(2*T))
    # The selected charged pion cap includes the final thermal factors.
    # For the forward bound use the min partner energy for threshold.
    pmax=np.sqrt((Emax-d)*(Emax+d))
    et=(4*mpi*mpi-2*d*d)/(2*(Emax+pmax))
    nba=min(n,hot*np.exp(-max(et,d)/(2*T)))
    up_pi=Kpi*nba*np.exp(-Emin/T)
    return np.array([up_low,up_high,up_pi])


def direct_mb_at_rest(d,T):
    """Independent on-shell Dirac spin trace scalar form, no final blocking."""
    def channel(mf):
        Emin=max(d,(4*mf*mf-2*d*d)/(2*d))
        pmin=np.sqrt((Emin-d)*(Emin+d))
        def f(p):
            E=np.hypot(d,p)
            s=2*d*(d+E)
            # v_Moller / beta_in, rationalized at p=0 using
            # (E-d)(E+d)=p^2. This exact rest formula is finite.
            vr=np.sqrt(s*(E+d)/(2*d))/E
            out=np.sqrt(max(0.,1-4*mf*mf/s))*(1+2*mf*mf/s)
            sigma_v=COUP2*s/(12*np.pi*((s-VECTOR_MASS**2)**2+
                        (VECTOR_MASS*VECTOR_GAMMA)**2))*(1+2*d*d/s)*out*vr
            return p*p*np.exp(-(E+d)/T)*sigma_v/np.pi**2
        if pmin>3.:
            # The extremely hard channel is irrelevant to this numerical
            # comparator. It remains covered in the analytic upper.
            return 0.
        return quad(f,pmin,3.,epsabs=1e-42,epsrel=2e-9)[0]
    return sum(channel(mm) for mm in LEPTONS)


def controls():
    errs=[]
    for d,T in ((m,.01),(M,.01),(M,.03)):
        n,mean=bath(d,T)
        rule=lambda i:quad(lambda p:p*p*np.hypot(d,p)**i*np.exp(-np.hypot(d,p)/T),
                           0.,3.,epsabs=1e-22,epsrel=1e-13)[0]/np.pi**2
        errs.extend([abs(rule(0)/n-1),abs(rule(1)/(n*mean)-1)])
    assert max(errs)<2e-13
    return max(errs)


def run(frozen,pion,stages):
    R,J=map(read,(frozen,pion))
    caps={m["q"]:max(d["pion_pair_sigma_v_cap_with_stated_support"] for d in m["controls"]) for m in J["cases"]}
    L=next(z for z in R["cases"] if z["q"]==32)
    assert len(stages)==len(L["stages"])
    rows=[]
    for j,S in enumerate(L["stages"]):
        u=S["fine"];T=S["T_GeV"]
        H=u["H_over_rest_decay_width"]*width
        p=u["cold_shell_edge_GeV"]; EE=np.hypot(p,m)
        lo,hi=u["local_scalar_energy_window_GeV"]
        low_H=m+lo; high_H=EE+hi
        v=upper(m,T,m,EE,caps[32])/H
        V=upper(M,T,low_H,high_H,caps[32])
        # A stationary source divided by the decay damping. With frozen
        # daughter distributions, damping >= width*M/E*(1-f_L,max).
        # Exact entropy reconstruction maximum in each archived shell.
        C=stages[j].get("checkpoint",stages[j])
        s=C["snapshot"]
        assert C["snapshot_sha256"]==S["snapshot_sha256"]
        a,b=np.asarray(s["entropy_variable_fermion"])
        fmax=float(expit(-min(np.min(a),np.min(a+b))))
        a,b=np.asarray(s["entropy_variable_scalar"])
        Bmax=float(1/np.expm1(min(np.min(a),np.min(a+b))))
        damp=width*M/high_H*(1-fmax)
        h=V/damp
        induced=(M/m)**3*width/H*(1+Bmax)*h
        f0=np.asarray(u["net_cold_shell_rate_per_Hubble_range"])
        extra=float(v.sum()+induced.sum())
        direct=direct_mb_at_rest(m,T)/H
        assert 0<direct<sum(v[:2])
        assert f0.max()+extra<0.
        rows.append(dict(T_GeV=T,base_sampled_net_range_per_Hubble=f0.tolist(),
                         direct_light_injection_upper_per_Hubble=v.tolist(),
                         independent_direct_light_lepton_MB_comparator_at_rest_over_Hubble=direct,
                         frozen_heavy_stationary_increase_upper=h.tolist(),
                         secondary_light_injection_upper_per_Hubble=induced.tolist(),
                         sampled_upper_net_including_injection_per_Hubble=float(f0.max()+extra),
                         max_added_to_smallest_sampled_base_sink=float(extra/(-f0.max())),
                         maximum_frozen_light_occupation=fmax))
    return dict(scope="Positive vacuum charged leptons plus supported charged-pion pair "
                      "injection (chosen cap); excludes unmodeled visible channels, elastic, "
                      "other internal reactions and all dynamic heavy transients. The true "
                      "light sign for the full coupled evolution is not inferred.",
                sequence_for_three_bounds=["leptons_below_sqrt_s_1","leptons_above_sqrt_s_1","selected_charged_pions"],
                independent_MB_quadrature_relative_defect=controls(),rows=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("frozen","pion","continued","observed"):p.add_argument(k,type=Path)
    a=p.parse_args()
    # Keep the same saved snapshots for the global instantaneous f_L peak.
    stages=read(a.continued)+read(a.observed)["blocks"]
    print(json.dumps(run(a.frozen,a.pion,stages),indent=2,allow_nan=False))
