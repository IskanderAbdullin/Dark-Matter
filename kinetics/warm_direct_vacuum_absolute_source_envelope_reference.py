"""Absolute scalar abundance headroom for the restricted direct vacuum sources."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import kve
from archived_collisionless_streaming_reference import Background
from conditional_relic_normalization_reference import conversion
from kloe_pion_band_portal_reference import table,NB_PER_GEV2
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import M,mp
from timelike_chemical_heat_eigenmodes_reference import leptons
from timelike_portal_onsager_coverage_reference import V,width,factor,gprod,mpi
from warm_babar_dispersive_gap_reference import POINT,mu_kernel_lower
from warm_late_cold_visible_offpole_dose_reference import log_tail
from warm_late_cold_visible_recoil_box_upper_reference import log_slope_minimum
from warm_light_occupation_stage_audit_reference import read


def clock(bg,n,tmax=.05,tmin=.0001,lower=False):
    assert 0 < tmin < tmax and bg.u[0] <= np.log(tmin) < np.log(tmax) <= bg.u[-1]
    T=np.geomspace(tmax,tmin,n+1)
    A=np.array([bg.scale(t) for t in T])
    t=T[:-1] if lower else T[1:]
    z=np.log(t)
    E=(2*np.pi**2/45)*t**5*np.exp(bg.gs(z)+.5*bg.gr(z))*np.sqrt(np.pi**2/90)/MPL
    return (T[1:] if lower else T[:-1]),np.diff(np.log(A))/E


def weight(lo,hi,T,U,bose,mass=M):
    assert np.all(lo>=4*mass*mass*(1-2e-14))
    den=(np.clip(V*V,lo,hi)-V*V)**2+(V*width)**2
    eu=np.sqrt(lo)[None,:]/T[:,None]
    G=T[:,None]*kve(1,eu)*np.exp(-eu)
    if bose:
        G*=(-1/np.expm1(-mpi/T))[:,None]**2
    p=2*factor*(hi+2*mass*mass)*hi**1.5*np.sqrt(1-4*mass*mass/hi)/den
    return (U@G)*p


def bounds(bg,tablefile,old,n,tmax=.05,tmin=.0001,mass=M):
    T,U=clock(bg,3*n,tmax,tmin)
    s=np.unique(np.r_[np.geomspace(4*mass*mass,16,n+1),[4*x*x for x in leptons if 2*mass<=2*x<=4]])
    l,h=s[:-1],s[1:]
    W=weight(l,h,T,U,False,mass)
    a=np.array(leptons)[:,None]**2/h
    R=(1+2*a)*np.sqrt(np.maximum(0,1-4*a))
    lepton=float((W*(h-l))@R.sum(axis=0))
    # A unit spectral transfer uses no more than an individual massless
    # charged-lepton weight; use the loose full three-lepton tail bound.
    pole=(np.hypot(V*V,V*width)+V*V)/(2*(V*width)**2)
    coef=2*12*(3*gprod*gprod*pole/(8*np.pi))/(120*np.pi**4)
    tail=np.exp(np.log(coef)+6*np.log(T)+np.array([log_tail(4/t) for t in T]))
    lepton+=float(U@tail)
    d,x,e=table(tablefile)
    few=5 if n<192 else 12
    sub=(d[:,None]-.005+.01*np.linspace(0,1,few+1))
    lo,hi=sub[:,:-1].ravel(),sub[:,1:].ravel()
    q=weight(lo,hi,T,U,True,mass)*hi/POINT
    kloe=float(np.sum(q*(hi-lo)*np.repeat(x+2*e,few)/NB_PER_GEV2))
    low=[]
    for row in old["bands"]:
        q=np.linspace(*row["s_GeV2"],n//4+1)
        lo,hi=q[:-1],q[1:]
        f=weight(lo,hi,T,U,True,mass)*hi/POINT
        H=mu_kernel_lower(hi)/(4*np.pi**3)
        low.append(float(max(f/H)*row["assumed_positive_dispersion_cap"]))
    q=np.geomspace(.95,16,n//2+1)
    b=float((np.diff(q)*weight(q[:-1],q[1:],T,U,True,mass)).sum()+U@(tail/(-np.expm1(-mpi/T))**2))
    return dict(invariant_base_bins=n,temperature_panels=len(T),lepton_yield_upper=lepton,
                KLOE_fixed_bin_plus_two_marginals_upper=kloe,additional_conditional_low_moment_upper_each=low,
                high_yield_upper_per_unit_R_pi=b,total_known_plus_unit_high_yield_upper=lepton+kloe+sum(low)+b)


def run(eos,kloe,gap,psd):
    b=Background(eos)
    slopes=[log_slope_minimum(b.gs,3),log_slope_minimum(b.gr,4)]
    prior=read(gap)
    baseline=float(sum(read(psd)["controls"][-1]["conditional_total_scalar_yield"]))
    results=[bounds(b,kloe,prior,n) for n in (96,240)]
    assert results[-1]["total_known_plus_unit_high_yield_upper"]>baseline
    assert results[-1]["total_known_plus_unit_high_yield_upper"]<results[0]["total_known_plus_unit_high_yield_upper"]
    v=conversion(eos)["Omega_h2_per_GeV_times_yield"]*mp
    C=results[-1]
    return dict(scope="Only the direct positive vacuum charged-lepton and conditional charged-pion sources, "
                      "50–.1 MeV adiabatic radiation. Fixed KLOE bin choices and moment caps are not confidence "
                      "limits; high pion contribution proportional to the separate assumed uniform R cap. "
                      "No other hadrons or nonlinear dark conversion/production bounded.",
               monotone_slopes=slopes,numerical_direct_lepton_and_KLOE_yield_12_to_50=baseline,Omega_h2_per_scalar_Y=v,
               direct_numerical_Omega_h2=baseline*v,direct_known_plus_unit_high_Omega_h2_upper=C["total_known_plus_unit_high_yield_upper"]*v,
               ratio_of_imposed_half_dark_omega_to_numerical_direct= .060/(baseline*v),
               rows=results)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","kloe","gap","psd"):
        p.add_argument(k,type=Path)
    a=p.parse_args();print(json.dumps(run(a.eos,a.kloe,a.gap,a.psd),indent=2,allow_nan=False))
