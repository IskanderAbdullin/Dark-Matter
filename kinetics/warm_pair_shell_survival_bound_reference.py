"""Shellwise lower survival from frozen parent profiles and bootstrapped Bose cap."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import expit

from scalar_pair_expanding_continuation_reference import moments
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import m
from warm_pair_selective_loss_bound_reference import loss_floor, moment_bounds, sha
from warm_scalar_pair_late_support_bound_reference import coefficient


def shell_bound(p,theta,n):
    K=p*p/(np.hypot(m,p)+m)
    x=np.linspace(0.,1.,n+1)
    k=K[:-1,None]+np.diff(K)[:,None]*x
    q=np.sqrt(k*(2*m+k))
    f=expit(-theta[0,:,None]-theta[1,:,None]*x)
    # Stable reconstruction aggregates spin/charge degeneracy g=4.
    # Use a stable difference-of-cubes expression for very fine panels.
    q0,q1=q[:,:-1],q[:,1:]
    vol=4/(6*np.pi**2)*(q1-q0)*(q1*q1+q1*q0+q0*q0)
    l=np.minimum(f[:,:-1],f[:,1:])*vol
    u=np.maximum(f[:,:-1],f[:,1:])*vol
    return [z.ravel() for z in (q1,l,u)]


def exposure(p,N,H,F):
    s=p*p
    results=[]
    k=(1.+F)**2/(1536*np.pi*m*m*H)
    for a in range(0,len(p),92):
        b=s[a:a+92,None]+s
        C=coefficient(np.sqrt(b/2.))
        results.extend((k*((C*b)@N)).tolist())
    return np.asarray(results)


def run(shape,audit,caps):
    s=json.loads(shape.read_text())["snapshots"][-1]
    d=json.loads(audit.read_text())["chained_endpoint_cases"][-1]
    v=json.loads(caps.read_text())
    assert s["T_visible_GeV"]==d["temperature_GeV"] and s["steps"]==d["n_steps"]
    z=next(x for x in v["cases"] if x["source_shape_sha256"]==sha(shape))
    T=s["T_visible_GeV"]
    assert z["anchor_T"]==T
    p=np.asarray(s["lab_momentum_edges_GeV"])
    th=[np.asarray(s[f]) for f in ("entropy_variable_scalar","entropy_variable_fermion")]
    N,U=moment_bounds(p,th)
    F=z["global_Bose_ceiling"]
    H=np.sqrt(np.pi**2/45)/MPL*(v["photon_floor_ratio"]*T)**2
    P=p[-1]
    assert P==z["finite_incoming_pmax_GeV"]
    # Old outside stable counts are not collision partners and are added once.
    E=d["all_tail_Omega_h2"][1]
    total=d["retained_and_all_particle_tail_Omega_h2"][1]
    kept=total-E
    assert 0.<E<total
    tau=N*coefficient(P)*P*P*(1+F)**2/(384*np.pi*m*m*H)
    r=U[-1]/(2*P*P)
    old=E+kept*loss_floor(tau,r)
    phases=[]
    for n in (4,16,64,256):
        highs,lo,hi=shell_bound(p,th[1],n)
        assert lo.sum()<=N<=hi.sum()
        x=exposure(highs,hi,H,F)
        assert np.min(x)>0 and np.max(x)<=tau*(hi.sum()/N)*(1+1e-12)
        # Normalize by the true total moment, not by lower or upper bins.
        obs=float(lo@np.exp(-x)/N)
        assert 0.<obs<=1.
        # Difference form bounds the loss much more effectively than
        # replacing the known total count by its lower quadrature bound.
        maxloss=float(-(hi@np.expm1(-x))/N)
        body=1.-maxloss
        assert body>=obs-1e-13
        phases.append(dict(kinetic_subpanels_per_shell=n,lower_and_upper_number_over_exact=[lo.sum()/N,hi.sum()/N],
                           maximum_exposure=float(np.max(x)),direct_lower_surviving_fraction=obs,
                           loss_upper_completion_lower_surviving_fraction=body,
                           inclusive_stable_Omega_floor=E+kept*body))
    assert np.all(np.diff([x["inclusive_stable_Omega_floor"] for x in phases])>=0.)
    return dict(input_sha256=[sha(a) for a in (shape,audit,caps)],photon_temperature_floor_ratio=v["photon_floor_ratio"],
        Bose_cap=F,source_active_stable_n=N,input_inclusive_stable_Omega=total,inert_old_tail_Omega=E,
        original_closed_Bose_moment_and_support_floor=old,panels=phases,
        scope="Equal unpolarized opposite charge occupations, ballistic loss-only parents on expanding finite support, no new parents, redistribution, inverse, inherited exterior scalars or re-entry. Newly produced scalar modes anywhere are allowed under the closed Bose cap. Not a bound on the fully coupled physical model.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for name in ("shape","audit","caps"):
        p.add_argument(name,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.shape,a.audit,a.caps),indent=2,allow_nan=False))
