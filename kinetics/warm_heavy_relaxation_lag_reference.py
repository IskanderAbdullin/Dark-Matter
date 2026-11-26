"""Exact heavy relaxation on a piecewise-linear saved-state surrogate."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import roots_legendre

from frw_exterior_restart_reference import import_restart
from joint_visible_heavy_adiabatic_reference import parent_branch,snapshot_functions
from late_radiation_background_reference import MPL,TabulatedRadiation
from three_species_decay_onsager_reference import M,width
from warm_isolated_recycling_energy_envelope_reference import moments
from warm_light_occupation_stage_audit_reference import read


def update(x,a,b,L):
    # Stable drift form of the driven relaxation equation over a cell.
    z=-np.expm1(-L)
    return b+(x-a)*np.exp(-L)-(b-a)*z/L


def independent():
    err=0.
    for L in (.03,.3,3.,30.):
        x,a,b=.012,.003,.043
        t=solve_ivp(lambda t,u: L*(a+(b-a)*t-u),[0.,1.],[x],
                    method="DOP853",rtol=2e-13,atol=1e-14).y[0,-1]
        err=max(err,abs(t-update(x,a,b,L)))
    assert err<1e-13
    return err


def run_case(eos,origin,tape,n):
    seed=read(origin)["snapshot"]
    stages=[seed]+[x["snapshot"] for x in read(tape)]
    A=seed["scale_factor"]
    B=np.array([s["scale_factor"]/A for s in stages])
    # Apply over a compact comoving band, fixed across the saved stages.
    s=np.linspace(0.,1.,101)
    x,w=roots_legendre(n)
    a,b=(s[:-1,None]+s[1:,None])/2,np.diff(s)[:,None]/2
    P=(a+b*x).ravel(); J=(b*w).ravel()
    F=[];rate=[];normal=[];boost=[]
    for R,stage in zip(B,stages,strict=True):
        T=stage["T_visible_GeV"]
        f,sp,g,fp=snapshot_functions(stage)
        p=P/R
        E=np.hypot(M,p)
        h=parent_branch(E,f,g,sp,fp,angle=13)
        total=moments(stage).sum(axis=0)[1]
        ext=import_restart(stage["positive_exterior_restart"]).energy(1.)
        H=np.sqrt((TabulatedRadiation(eos,T).state(1.)[1]+total+ext)/3)/MPL
        lam=width*M*h[:,-1]/(E*H)
        F.append(h[:,0]);rate.append(lam);normal.append(total);boost.append(E)
    y=F[0].copy()
    rows=[]
    for j in range(1,len(B)):
        lam=np.minimum(rate[j-1],rate[j])
        L=lam*np.log(B[j]/B[j-1])
        y=update(y,F[j-1],F[j],L)
        d=y-F[j]
        om=lambda z: float((2/np.pi**2)*np.sum(z*boost[j]*P**2*J)/normal[j]/B[j]**3)
        rows.append(dict(b_relative=float(B[j]),min_integrated_relaxation=float(min(L)),
                         peak_abs_lag=float(max(abs(d))),energy_abs_lag_over_core=om(abs(d)),
                         stationary_E_fraction=om(F[j])))
    return rows


def run(eos,p16,p32,t16,t32):
    rows=[]
    for q,p,t in ((16,p16,t16),(32,p32,t32)):
        a,b=[run_case(eos,p,t,n) for n in (4,10)]
        for x,y in zip(a,b,strict=True):
            for k in ("energy_abs_lag_over_core","stationary_E_fraction"):
                assert abs(x[k]/y[k]-1.)<2e-3
        rows.append(dict(q=q,checks=[a,b]))
    return dict(scope="Thinned surrogate: between each saved pair, linearly interpolate "
                     "the stationary source in log expansion and freeze its relaxation "
                      "coefficient at the smaller endpoint value. Starts at stationary "
                      "heavy population and omits extra heavy sources. Exact for that "
                      "interpolation, not an error bound on actual intermediate shapes, "
                      "unrecorded exterior, or the complete heavy evolution.",
                elementary_update_absolute_test_error=independent(),cases=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("eos","p16","p32","t16","t32"):
        p.add_argument(name)
    a=p.parse_args()
    print(json.dumps(run(*vars(a).values()),indent=2,allow_nan=False))
