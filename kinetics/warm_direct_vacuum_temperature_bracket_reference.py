"""Positive direct candidate abundance enclosures as prescribed reheating changes."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.special import kve, roots_legendre

from empty_lepton_pair_spectrum_reference import Source
from sm_dark_entropy_exchange_reference import MPL
from warm_direct_vacuum_scalar_psd_reference import pion_parents

from warm_direct_vacuum_absolute_source_envelope_reference import (
    Background, M, V, width, factor, leptons, table, NB_PER_GEV2, POINT,
    bounds, clock, log_slope_minimum, read)


def geometry(l,h,T):
    """Lower raw MB matrix over invariant panels (without bath factors)."""
    den=np.maximum((l-V*V)**2,(h-V*V)**2)+(V*width)**2
    beta=np.sqrt(np.maximum(0,(l-4*M*M)/l))
    x=np.sqrt(h)[None,:]/T[:,None]
    return (2*factor*(l+2*M*M)*l**1.5*beta/den)[None,:]*T[:,None]*kve(1,x)*np.exp(-x)


def lower(bg,kloe,n,tmax,tmin):
    tc,U=clock(bg,3*n,tmax,tmin,lower=True)
    th=np.geomspace(tmax,tmin,3*n+1)[:-1]
    s=np.unique(np.r_[np.geomspace(4*M*M,16,n+1),[4*m*m for m in leptons if .2<=2*m<=4]])
    l,h=s[:-1],s[1:]
    A=geometry(l,h,tc)
    pairs=[]
    for m in leptons:
        a=m*m/l
        R=(1+2*a)*np.sqrt(np.maximum(0,1-4*a))
        # E1+E2>=sqrt(s), Ei>=m; maximum FD denominator occurs at
        # maximal energy imbalance and its smallest admissible sum.
        stat=1/((1+np.exp(-m/th[:,None]))*(1+np.exp(-(np.maximum(np.sqrt(l)[None,:]-m,m))/th[:,None])))
        # Below threshold, this channel is exactly zero. Using the
        # clipped kinematic expression avoids negative square roots.
        pairs.append(float(U@((A*stat)@(R*(h-l)))))
    d,x,e=table(kloe)
    sub=d[:,None]-.005+.01*np.linspace(0,1,13)
    l,h=sub[:,:-1].ravel(),sub[:,1:].ravel()
    # BE >= MB on the incoming charged-pion legs. No low moment
    # or high continuum is used in this lower comparison.
    pion=float(U@(geometry(l,h,tc)@((h-l)*l/POINT*np.repeat(np.maximum(0,x-2*e),12)/NB_PER_GEV2)))
    return pairs,pion


def numerical_central(bg,kloe,edges,order):
    """Independent lab-energy/angle fold, not an endpoint enclosure."""
    data=table(kloe)
    src=Source(((9,28,24),(14,46,32))[order],tail=75.,dark_mass=M)
    x,w=roots_legendre((20,34)[order])
    counts=np.zeros((len(edges)-1,2))
    for i,(l,r) in enumerate(zip(np.log(edges[:-1]),np.log(edges[1:]),strict=True)):
        for u,du in zip((r+l)/2+(r-l)*x/2,(r-l)*w/2,strict=True):
            t=np.exp(u)
            sh=(2*np.pi**2/45)*t**5*np.exp(bg.gs(u)+.5*bg.gr(u))*np.sqrt(np.pi**2/90)/MPL
            clock_factor=du*(1+bg.gs.derivative()(u)/3)/sh
            lep=2*src.event_moments(t)[3].sum()
            pi=pion_parents(t,data,*((5,24,18),(8,36,26))[order])[1].sum()
            counts[i]+=clock_factor*np.asarray((lep,pi))
    return counts.cumsum(axis=0)


def run(eos,kloe,dispersion,reference):
    bg=Background(eos)
    v=read(reference)
    prior=read(dispersion)
    unit=v["Omega_h2_per_scalar_Y"]
    slope=[log_slope_minimum(bg.gs,3),log_slope_minimum(bg.gr,4)]
    ladder=(.025,.04,.045,.05,.0555,.06,.075,.1)
    recent=tuple(.9*t for t in ladder)
    edges=np.array(sorted({.006,.012,*ladder,*recent}))
    nc=[numerical_central(bg,kloe,edges,j)*unit for j in (0,1)]
    rows=[]
    for t in ladder:
        controls=[]
        for n in (120,240,480):
            c=bounds(bg,kloe,prior,n,t)
            a,p=lower(bg,kloe,n,t,.0001)
            le=float(sum(a))
            lo=(le+p)*unit
            up=c["total_known_plus_unit_high_yield_upper"]*unit
            assert 0<lo<up
            controls.append(dict(invariant_panels=n,lower=lo,upper=up,
                lepton_lower_by_species=np.asarray(a).tolist(),KLOE_lower=p,
                high_pion_Omega_per_unit_R_upper=c["high_yield_upper_per_unit_R_pi"]*unit,
                upper_terms=c))
        assert all(x["lower"]<=y["lower"] and x["upper"]>=y["upper"] for x,y in zip(controls,controls[1:]))
        if t==.05:
            assert abs(controls[1]["upper"]/v["direct_known_plus_unit_high_Omega_h2_upper"]-1)<5e-13
            assert controls[-1]["lower"]<v["direct_numerical_Omega_h2"]<controls[-1]["upper"]
        # Compare the positive recent dose with the independent upper
        # for the same modeled channels over the whole history.
        window=.9*t
        wl,wp=lower(bg,kloe,480,t,window)
        recent_lower=(sum(wl)+wp)*unit
        old=bounds(bg,kloe,prior,480,window)
        old_upper=old["total_known_plus_unit_high_yield_upper"]*unit
        j=np.searchsorted(edges,t)-1
        h=np.searchsorted(edges,window)-1
        total=float(nc[1][j].sum())
        frac=float((nc[1][j]-nc[1][h]).sum()/total)
        assert controls[-1]["lower"]<total<controls[-1]["upper"]
        assert max(abs(nc[1][j]/nc[0][j]-1.))<1e-3
        if t==.05:
            j12=np.searchsorted(edges,.012)-1
            assert abs((nc[1][j]-nc[1][j12]).sum()/v["direct_numerical_Omega_h2"]-1.)<3e-4
        rows.append(dict(T_reheat_GeV=t,controls=controls,
                         central_lepton_and_KLOE_Omega_h2=nc[1][j].tolist(),
                         central_lower_fold_cutoff_GeV=edges[0],
                         central_max_relative_quadrature_change=float(max(abs(nc[1][j]/nc[0][j]-1.))),
                         central_fraction_from_top_tenth=frac,
                         dose_above_T_GeV=window,
                         fraction_attributed_to_top_tenth_lower=min(1.,max(
                            recent_lower/controls[-1]["upper"],
                            1-old_upper/controls[-1]["lower"]))))
    assert all(x["controls"][-1]["lower"]<y["controls"][-1]["lower"]
               and x["controls"][-1]["upper"]<y["controls"][-1]["upper"] for x,y in zip(rows,rows[1:]))
    return dict(scope="Exact positive panel enclosures conditional on instantaneous specified reheating, "
                      "adiabatic visible radiation, vacuum direct thermal lepton and charged-pion channels, "
                      "two scalar daughters per heavy-pair event, and free survival. Entire lower omits unmeasured "
                      "pion windows. Pion upper uses two given positive moments and the assumed high-band "
                      "R_pi<=1; measured bin choices +/- two stated marginal errors are not a statistical interval. "
                      "No stable prediction, nonlinear recycling, elastic heating, or other number-changing source.",
                eos_sha256=hashlib.sha256(eos.read_bytes()).hexdigest(),
                KLOE_sha256=hashlib.sha256(kloe.read_bytes()).hexdigest(),
                log_slope_lower= slope, Y_to_Omega_h2=unit, rows=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("eos","kloe","dispersion","reference"):
        p.add_argument(key,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.kloe,a.dispersion,a.reference),indent=2,allow_nan=False))
