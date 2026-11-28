"""Tighten the direct-lepton envelope by separating sub-pole total energies."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.special import gammaln,logsumexp
from archived_collisionless_streaming_reference import Background
from decay_scalar_ir_relaxation_reference import limits
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import M,P
from timelike_portal_onsager_coverage_reference import V,width,gprod
from warm_late_cold_visible_vacuum_dose_reference import clock_monotonicity,run as base


def log_tail(x):
    j=np.arange(6)
    return gammaln(6)-x+logsumexp(j*np.log(x)-gammaln(j+1))


def log_energy_bin(lo,hi,T):
    a,b=log_tail(lo/T),log_tail(hi/T) if np.isfinite(hi) else -np.inf
    assert a>b
    return a+np.log(-np.expm1(b-a))


def source(T,q,edges):
    xmin=limits(q)[1]+M
    b=V*width
    pole=np.hypot(V*V,b)
    a=[]
    for lo,hi in zip(edges[:-1],edges[1:],strict=True):
        floor=max(lo,xmin)
        if hi<=floor:
            continue
        # s cannot exceed total energy squared. The BW shape s / BW
        # increases up to sqrt(V^4+b^2), then decreases. Ignoring a
        # lower invariant-mass bound is conservative.
        s=min(pole,hi*hi)
        k=s/((s-V*V)**2+b*b)
        a.append(np.log(k)+log_energy_bin(floor,hi,T))
    coeff=12*(3*gprod*gprod/(8*np.pi))/(120*np.pi**4)
    return np.log(coeff)+6*np.log(T)+logsumexp(a)


def run(eos,archive):
    old=base(eos,archive)
    bg=Background(eos)
    clock_monotonicity(bg)
    y=old["comparison_archive_held_equality_cold_yield"]
    low=np.array([0,.5,1,2,4,7,9,10,np.inf])
    fine=np.unique(np.r_[low,np.linspace(0,10,81)])
    out=[]
    for a in old["rows"]:
        hi,lo=a["T_max_GeV"],a["T_min_GeV"]
        T=np.geomspace(hi,lo,513)
        scale=np.array([bg.scale(t) for t in T])
        p=.1*bg.T0/scale
        time=np.log(2*np.diff(np.log(scale)))
        denom=np.log(np.sqrt(np.pi**2/45)*T[1:]**2/MPL*4*np.pi**2/45*T[1:]**3)
        angle=np.log(p[:-1]**2/(2*P*(P+np.sqrt((P-p[:-1])*(P+p[:-1])))))
        control=[]
        for edges in (low,fine):
            src=[source(t,q,edges) for t,q in zip(T[:-1],p[:-1],strict=True)]
            integrand=time+src-denom
            Y=float(logsumexp(integrand))
            YA=float(logsumexp(integrand+angle))
            control.append([float((z-np.log(y))/np.log(10.)) for z in (Y,YA)])
        ref=a["refinements"][-1]["log10_upper_over_archive_cold_yield"]
        assert control[1][0]<=control[0][0]<=ref+1e-12
        assert control[1][1]<=control[0][1]<=control[0][0]
        out.append(dict(T_start_GeV=hi,log10_previous_over_held_cold=ref,
                    log10_total_energy_bin_bounds_over_held_cold=control))
    return dict(columns=["log10_kinematic_upper_over_hold","log10_with_open_unpolarized_decay_angle_over_hold"],scope="Same vacuum direct thermal charged-lepton envelope; subdivide total incoming energies "
                     "to avoid granting the far-away resonance maximum to all thermal collisions.",rows=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("archive",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.archive),indent=2,allow_nan=False))
