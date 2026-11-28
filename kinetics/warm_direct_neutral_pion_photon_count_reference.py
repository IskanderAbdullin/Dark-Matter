"""Transfer a positive radiative neutral-pion moment, treating the soft photon explicitly."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import kve

from archived_collisionless_streaming_reference import Background
from three_species_decay_onsager_reference import M,m
from timelike_portal_onsager_coverage_reference import alpha,V,width,factor
from warm_direct_vacuum_absolute_source_envelope_reference import clock,POINT,mu_kernel_lower,read
from warm_neutral_pion_photon_budget_reference import m0

# Conditional full positive moment cap from the published data-constrained dispersion fit.
# It can be spent on a strict subset; no actual low-energy rate is inferred from doing so.
CAP=(43.8+2*.6)*1e-11
# Independent explicit factor-two allowance above the fit's sub-.60 estimate.
LOW_CAP=2*1.3e-11


def pieces(l,h,T,U):
    """Clock-folded upper MB and one-energy-weight MB terms before multiplying R."""
    q=T[:,None]
    r=np.sqrt(l)[None,:]
    d=np.sqrt(h)[None,:]
    x=r/q
    k0,k1=(kve(i,x)*np.exp(-x) for i in (0,1))
    w=U/(-np.expm1(-m0/T))
    n0=w@(q*k1)
    n1=w@(4*q*(2*q*q*k1+q*d*k0))
    b=(np.clip(V*V,l,h)-V*V)**2+(V*width)**2
    rho=lambda mass:(h+2*mass*mass)*np.sqrt(np.maximum(0.,1-4*mass*mass/h))
    K=2*factor*h**1.5*(2*rho(M)+rho(m))/b
    return n0*K,n1*K


def low_window(split,T,U,n):
    s=np.linspace(m0*m0,split*split,n+1)
    l,h=s[:-1],s[1:]
    a,b=pieces(l,h,T,U)
    d=np.maximum(h-m0*m0,0.)
    # Here R proportional to d^3/s^2 eliminates the massless singularity
    # analytically, so do not estimate 1/(s-m0²) on the first panel.
    r=np.pi*alpha/2*d**3/h**2
    r_div=np.pi*alpha/2*(d/h)**2
    return float(np.diff(s)@(r*a+r_div*b))


def distant(a,b,T,U,n):
    s=np.geomspace(a*a,b*b,n+1)
    l,h=s[:-1],s[1:]
    x,y=pieces(l,h,T,U)
    return s,(x+y/(l-m0*m0))


def run(eos,ladder,base):
    units=read(ladder)["Y_to_Omega_h2"]
    B=read(base)
    bg=Background(eos)
    out=[]
    for row in B["rows"]:
        T0=row["T_reheat_GeV"]
        prior=max(t["total_upper_including_only_measured_bands"] for t in row["controls"][-1]["control_cases"])
        tiers=[]
        for n in (240,480):
            T,U=clock(bg,3*n,T0)
            sid,high=distant(1.35,4.,T,U,n)
            above=units*float(np.diff(sid)@high)
            cross=[]
            for v in (.32,.42,.60):
                s,W=distant(v,1.35,T,U,n)
                den=mu_kernel_lower(s[1:])/(4*np.pi**3)
                ratios=W*s[1:]/(POINT*den)
                imax=int(np.argmax(ratios))
                upper=units*CAP*float(ratios[imax])
                lower=units*low_window(v,T,U,n)
                cross.append(dict(separation_GeV=v,low_omega_ceiling_per_GeVminus2_Fhat_squared=lower,
                    upper_count_omega_from_positive_moment=upper,invariant_controlling_panel_GeV2=s[imax:imax+2].tolist(),
                    total_omega_if_low_cap_Fhat1_and_no_others=prior+upper+lower,
                    total_omega_if_low_cap_Fhat10_and_no_others=prior+upper+100*lower))
            low_coefficient=cross[0]["low_omega_ceiling_per_GeVminus2_Fhat_squared"]
            s,w=distant(.32,.60,T,U,n)
            sub=units*LOW_CAP*float(max(w*s[1:]/(POINT*mu_kernel_lower(s[1:])/(4*np.pi**3))))
            high=cross[-1]["upper_count_omega_from_positive_moment"]
            tied=dict(low_Fhat_squared_coefficient=low_coefficient,
                      mid_032_060_count_omega_upper=sub,high_060_135_count_omega_upper=high,
                      total_omega_with_Fhat10_local_and_disjoint_moments=prior+100*low_coefficient+sub+high)
            assert all(z["total_omega_if_low_cap_Fhat1_and_no_others"]>.0 for z in cross)
            tiers.append(dict(panels=n,cross=cross,disjoint_moment_variant=tied,tail_1p35_to_4_omega_per_unit_R=above))
        for a,b in zip(tiers[0]["cross"],tiers[1]["cross"],strict=True):
            assert b["upper_count_omega_from_positive_moment"]<a["upper_count_omega_from_positive_moment"]
            assert b["low_omega_ceiling_per_GeVminus2_Fhat_squared"]<a["low_omega_ceiling_per_GeVminus2_Fhat_squared"]
        out.append(dict(T_reheat_GeV=T0,prior_upper=prior,controls=tiers))
    return dict(scope="Exclusive pi0 photon vacuum slot only. Conditional positive moment through 1.35 GeV "
                      "taken from a data-constrained dispersive fit and spent on the upper interval alone. "
                      "Below each tested split assume a separate constant |Fhat| cap; overlapping split alternatives "
                      "are never summed. Massless photon enhancement bounded at fixed total lab energy Q. "
                      "No new channels above 1.35 included in quoted totals; a separate finite weight only "
                      "covers 1.35–4 GeV per unit R. No thermal photon masses, in-medium spectral states or "
                      "dark regeneration bound.",
                source="https://arxiv.org/abs/2007.12696",positive_a_mu_cap=CAP,sub_060_conditional_a_mu_cap=LOW_CAP,rows=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("eos","ladder","baseline"):
        p.add_argument(key,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.ladder,a.baseline),indent=2,allow_nan=False))
