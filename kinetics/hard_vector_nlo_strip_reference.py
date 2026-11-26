"""Strict perturbative vector trace at fixed invariant-mass strips only."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

from born_thermal_quark_hard_portal_reference import (
    UEDGES,gauss,kernel,end,thermal_trace_over_vacuum)


def spectra(x,k):
    e=np.hypot(x,k)
    u=k*k/(e+x)
    rho0=3*x*x/(4*np.pi)
    return u,rho0*thermal_trace_over_vacuum(x,u)


def strip(a,method,drop=False):
    x,k,Lo,c=a.T
    assert np.all(x==x[0]) and len(a)==len(UEDGES)
    u,L=spectra(x,k)
    assert np.max(abs(u-np.asarray(UEDGES)))<5e-13
    assert np.max(abs(Lo/L-1.))<5e-12
    j=np.arange(len(a))
    if drop:
        j=j[(j%2==0)|(j==len(a)-1)]
    X=x[0]
    if method=="pchip":
        interp=PchipInterpolator(u[j],c[j])
    elif method=="linear":
        def interp(t):
            return np.interp(t,u[j],c[j])
    else:
        raise ValueError(method)
    def fold(edges):
        z,w=gauss(edges,8)
        q=np.sqrt(z*(2*X+z))
        Lo2=spectra(X,q)[1]
        d=interp(z)
        y=(q/X)/np.expm1(X+z)*w
        moments=np.array([y,y*(X+z)])
        assert np.all(d>0.)
        return moments@Lo2,moments@d
    B,G=fold(u)
    ratios=G/B
    # Optional lab-frame buffer: (omega-k)/T>=1 restricts
    # u to (x-1)^2/2 at the fixed invariant mass.
    U=(X-1.)**2/2
    mask=(u[0]<U<u[-1])
    assert mask
    kept_B,kept_G=fold(np.unique(np.minimum(u,U)))
    return dict(x_invariant_mass_over_T=float(X),u_boost_bounds=[float(u[0]),float(u[-1])],
                interpolation=method,probes_used=len(j),
                alpha_s_variations=[dict(alpha_s=al,strict_NLO_part_over_thermal_Born_rate_and_heat=
                                  (4*np.pi*al*ratios).tolist()) for al in (.2,.3,.4)],
                additional_omega_minus_k_over_T_ge_1_cut=dict(
                     u_max=float(U),kept_Born_count_and_heat_fraction=(kept_B/B).tolist(),
                       strict_NLO_part_over_thermal_Born_rate_and_heat_at_alpha0p3=
                                  (4*np.pi*.3*kept_G/kept_B).tolist()))


def report(data,focused,strict):
    r=np.loadtxt(data)
    f=np.atleast_2d(np.loadtxt(focused))
    u,c=spectra(f[:,0],f[:,1])
    err=float(np.max(abs(f[:,2]/c-1.)))
    assert err<5e-12
    check=[]
    if strict:
        s=np.atleast_2d(np.loadtxt(strict))
        for x,k,v,g in s:
            z=r[(r[:,0]==x)&(np.abs(r[:,1]-k)<1e-12)]
            assert len(z)==1
            check.append(dict(x=float(x),k_over_T=float(k),abs_native_Lo_difference=float(abs(v/z[0,2]-1.)),
                              abs_relative_native_strict_NLO_change=float(abs(g/z[0,3]-1.))))
    return dict(focused_probe_columns=["mass_over_T","k_over_T","rhoV_Born_over_T2","rhoV_NLO_coefficient_over_T2"],
                focused_probe_values=f.tolist(),
                focused_Born_trace_check_relative_defect=err,
                focused_strict_NLO_over_LO_at_alpha0p3=(4*np.pi*.3*f[:,3]/f[:,2]).tolist(),
                direct_tolerance_rechecks=check,
                complete_mass_strip_controls=[
                   strip(r[r[:,0]==x],mode,drop)
                   for x in (4.,4.5) for mode,drop in (("pchip",False),("linear",False),("pchip",True))],
                scope="Fixed selected invariant masses only. Integrate boost u, not invariant mass or cosmic temperature. Strict two-loop part for dimensionless fixed alpha_s, no LPM, lattice matching, running-coupling uncertainty, cosmological relic result or full spectral coverage.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("strips",type=Path)
    p.add_argument("focus",type=Path)
    p.add_argument("--strict",type=Path)
    p.add_argument("--eos",type=Path)
    a=p.parse_args()
    r=report(a.strips,a.focus,a.strict)
    if a.eos:
        from conditional_relic_normalization_reference import conversion
        from late_radiation_background_reference import TabulatedRadiation
        bg=TabulatedRadiation(a.eos,.05)
        f=conversion(a.eos)["Omega_h2_per_GeV_times_yield"]
        full=kernel(8)
        narrow=kernel(8,[4.,4.5])
        ks=("count_omega_if_all_convert","energy_floor_over_radiation_at_1MeV")
        summary=[]
        for lo,hi in ((.25,.3486245429),(.30,.42)):
            B=end(bg,f,full,lo,hi,8)
            N=end(bg,f,narrow,lo,hi,8)
            summary.append(dict(T_range_GeV=[lo,hi],
               pure_Born_FD_fraction_for_x_4_to_4p5_only=[N[k][0]/B[k][0] for k in ks]))
        r["coverage_not_correction"]=dict(fd_count_then_heat_fractions=summary,
            scope="Only pure Born coverage over the narrow mass range. Does not interpolate or propagate NLO coefficients.")
    print(json.dumps(r,indent=2,allow_nan=False))
