"""Hypothetical common two-bath equilibrium from archived total number and energy."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.special import expit,roots_legendre

from frw_exterior_restart_reference import import_restart
from late_radiation_background_reference import MPL,TabulatedRadiation
from three_species_decay_onsager_reference import m,mp
from warm_common_pair_onsager_reference import response
from warm_common_pair_thermal_relaxation_reference import modes
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom


def solve(target,cap,n):
    x,w=roots_legendre(n)
    E=np.r_[0.,np.geomspace(1e-7,cap,30)]
    h=np.diff(E)/2
    p=(E[:-1,None]+E[1:,None])/2+h[:,None]*x
    weight=(h[:,None]*w*p*p/(2*np.pi**2)).ravel()
    es=np.array([np.hypot(p.ravel(),z) for z in (mp,m)])
    G=np.array([1,4])[:,None]

    def eval_at(u):
        T=np.exp(u[0]); eta=u[1]
        q=es/T-eta
        f=np.array([1/np.expm1(q[0]),expit(-q[1])])
        W=G*f*weight
        return np.array([W.sum(axis=1),(W*es).sum(axis=1)])

    r=least_squares(lambda u:np.log(eval_at(u).sum(axis=1)/target),[np.log(.007),-7.5],
                    bounds=([np.log(.001),-100.],[np.log(.06),-1e-4]),
                    xtol=3e-14,ftol=3e-14,gtol=3e-14)
    T=np.exp(r.x[0]); mu=T*r.x[1]
    Z=eval_at(r.x)
    err=max(abs(Z.sum(axis=1)/target-1.))
    assert err<2e-12
    return float(T),float(mu),Z,float(err)


def report(eos,path,relax):
    record=read(path)
    record=record[-1] if isinstance(record,list) else record
    S=record["snapshot"];assert digest(S)==record["snapshot_sha256"]
    p=np.asarray(S["lab_momentum_edges_GeV"])
    N,K,E=get_mom(p,[S["entropy_variable_"+b] for b in ("scalar","fermion")])
    old=np.array([N.sum(axis=1),E.sum(axis=1)])
    target=old.sum(axis=1)
    out=[]
    for cap in (p[-1],1.07):
        T,mu,z,d=solve(target,cap,16)
        T2,M2,check,d2=solve(target,cap,24)
        np.testing.assert_allclose([T,mu],[T2,M2],rtol=1e-11)
        a=dict(cap_GeV=cap,common_temperature_GeV=T,common_mu_GeV=mu,max_total_relative_defect=max(d,d2),
               equilibrium_number_and_energy_by_species=z.tolist(),
               equilibrium_to_archived_light_number_energy=(z[:,1]/old[:,1]).tolist())
        if relax and cap==1.07:
            q=np.r_[0.,np.geomspace(.00004,.0015,6),np.geomspace(.003,cap,11)]
            raw=response(T,mu,q,(3,7,6,2),steps=(4e-5,1.2e-5))
            report=modes(raw,eos,18)
            visible=TabulatedRadiation(eos,S["T_visible_GeV"]).state(1.)[1]
            H=np.sqrt((visible+target[1]+import_restart(S["positive_exterior_restart"]).energy(1.))/3)/MPL
            a["synthetic_pair_tangent_rates_over_actual_anchor_H"]=(np.asarray(report["decay_rates_GeV"])/H).tolist()
            a["equilibrium_FD_BE_response"]=raw
        out.append(a)
    return dict(snapshot_sha256=digest(S),T_GeV=S["T_visible_GeV"],archived_species_number_energy=old.tolist(),controls=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("paths",nargs="+",type=Path)
    p.add_argument("--relax",action="store_true")
    a=p.parse_args()
    print(json.dumps(dict(scope="Hypothetical internally thermal reservoirs matching only the archived aggregate conserved moments.",
                          cases=[report(a.eos,p,a.relax) for p in a.paths]),indent=2,allow_nan=False))
