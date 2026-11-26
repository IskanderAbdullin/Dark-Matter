"""Uniform frozen tagged-F=1/2 loss margin for the scalar-pair term."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import expit

from three_species_decay_onsager_reference import M,m,mp,y
from warm_light_occupation_stage_audit_reference import digest,read
from warm_pair_cap_affine_event_reference import from_archive


def scalar_minima(E,Z,threshold):
    """Minimize one entropy and the sum over boxes intersected by E1+E2>=threshold."""
    lo=np.maximum(E[:-1],threshold-E[-1]); hi=E[1:]
    slope=Z[1]/np.diff(E)
    def val(j,x):return Z[0,j]+slope[j]*(x-E[j])
    keep=lo<=hi
    j=np.flatnonzero(keep)
    assert len(j)>0
    a0=np.min(np.r_[val(j,lo[j]),val(j,hi[j])])
    best=np.inf
    for j in np.flatnonzero(keep):
        for k in np.flatnonzero(keep):
            V=[(x,t) for x in (lo[j],hi[j]) for t in (lo[k],hi[k])]
            V += [(x,threshold-x) for x in (lo[j],hi[j])]
            V += [(threshold-t,t) for t in (lo[k],hi[k])]
            for x,t in V:
                if x+t>=threshold-1e-14 and lo[j]-1e-14<=x<=hi[j]+1e-14 and lo[k]-1e-14<=t<=hi[k]+1e-14:
                    best=min(best,val(j,x)+val(k,t))
    assert a0>0. and np.isfinite(best)
    return a0,best


def bound(path):
    data=from_archive(path)
    row=read(path); row=row[-1] if isinstance(row,list) else row
    S=row["snapshot"]
    assert digest(S)==data["snapshot_sha256"]
    p=np.asarray(S["lab_momentum_edges_GeV"])
    scalar=np.asarray(S["entropy_variable_scalar"])
    light=np.asarray(S["entropy_variable_fermion"])
    B,F=np.hypot(mp,p),np.hypot(m,p)
    et=F[-1]
    gap=data["min_energy_only_entropy_gap"]["gap"]
    assert gap>0.
    # With the tagged value fixed at one half, G/L=exp(-gap) per
    # event whenever the partner is represented. If any scalar is
    # exterior its vacuum occupation simply makes G=0.
    inside_factor=-np.expm1(-gap)
    s0=2*m*m+2*m*F[:-1]
    s1=2*m*m+2*(et*F[1:]+p[-1]*p[1:])
    ap=np.sqrt(np.maximum(0.,s0/4-m*m))
    As=M*M-m*m-mp*mp+s1/2
    # For incoming cos<=0, longitudinal term magnitude >=2 p* M/A.
    amin=8*y**4*ap**2*M*M/As**2
    flux=.5*np.sqrt(1-4*mp*mp/s0)
    fmin=expit(-np.maximum(light[0],light[0]+light[1]))
    per_bin=amin*flux*fmin*p[:-1]*np.diff(F)/(128*np.pi**3*et)
    # All these half-angular rectangles are disjoint partner shells.
    guaranteed_loss=float(per_bin.sum())*inside_factor
    # Exterior partner must satisfy E>=E_F,max, so both scalars must
    # lie in E_i>=m+E_F,max-E_scalar,max with sum at least m+E_F,max.
    # The vacuum partner can add gain but has no forward loss.
    emin,emax=et,2*B[-1]-m
    if emax<=emin:
        zmin=zsum=None
        gain=0.
    else:
        zmin,zsum=scalar_minima(B,scalar,m+emin)
        boson_product=np.exp(-zsum)/(-np.expm1(-zmin))**2
        upper_s=(2*B[-1])**2
        ps=np.sqrt(upper_s/4-m*m)
        C=M*M-m*m-mp*mp
        max_amp=2*y**4*((upper_s/4)*(2*ps/C)**2+(ps*(4*m+2*M)/C)**2)
        phase=(emax-emin)*np.sqrt((emax-m)*(emax+m))
        gain=float(max_amp*boson_product*phase/(128*np.pi**3*m))
    ratio=gain/guaranteed_loss
    assert ratio<1e-4
    return dict(scope="Uniform in supported tagged momentum at one frozen bath; vacuum partners outside. Scalar pair only.",
                T_GeV=S["T_visible_GeV"], snapshot_sha256=digest(S),
                entropy_gap_inside=gap,inside_maximum_gain_over_loss=float(np.exp(-gap)),
                exterior_scalar_minimum_entropy=zmin,exterior_scalar_pair_minimum_entropy_sum=zsum,
                interior_loss_lower_GeV=guaranteed_loss,exterior_gain_upper_GeV=gain,
                exterior_upper_over_interior_lower=ratio,net_upper_GeV=gain-guaranteed_loss)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("paths",nargs="+",type=Path)
    print(json.dumps([bound(f) for f in p.parse_args().paths],indent=2,allow_nan=False))
