"""Upper bound on bare high inverse decay from finite q64 shell moments alone."""

import argparse
import json
from pathlib import Path

import numpy as np

from factorized_prompt_return_graph_reference import C
from frw_exterior_restart_reference import import_restart
from late_radiation_background_reference import MPL,TabulatedRadiation
from three_species_decay_onsager_reference import M,m,mp
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_high_parent_ledger_reference import high
from warm_prompt_reentry_shape_audit_reference import get_mom


def one(s,p_cut,eos):
    p=np.asarray(s["lab_momentum_edges_GeV"])
    n,k,e=get_mom(p,[s["entropy_variable_scalar"],s["entropy_variable_fermion"]])
    mm=np.asarray([mp,m])
    E=np.hypot(p,mm[:,None])
    target=np.hypot(M,p_cut)
    H=np.sqrt((TabulatedRadiation(eos,s["T_visible_GeV"]).state(1.)[1]+e.sum()+
               import_restart(s["positive_exterior_restart"]).energy(1.))/3)/MPL
    N=n.shape[1]
    H0=target-E[1,1:][None,:]
    H1=target-E[0,1:][:,None]
    tests=[]
    for use_energy in (False,True):
        legs=[]
        for j,hh in enumerate((H0,H1)):
            axis=(slice(None),None) if j==0 else (None,slice(None))
            low=E[j,:-1][axis]; hi=E[j,1:][axis]; t=np.maximum(low,hh)
            shell_n=n[j][axis]
            if use_energy:
                excess=np.maximum(0.,(k[j]*(1+2e-12)-n[j]*(E[j,:-1]-mm[j]))[axis])
                frac=np.divide(excess,np.maximum(t-low,1e-200),out=np.zeros((N,N)),where=t>low)
                eligible=np.where(t>low,np.minimum(shell_n,frac),shell_n)
            else:
                eligible=shell_n
            t=np.where(t<=hi,t,hi)
            p0=np.sqrt(np.maximum(0.,(t-mm[j])*(t+mm[j])))
            denom=p0*t
            allowed=(hh<=hi)&(denom>0)
            # A shell pair with nonzero weights and adequate total energy
            # always puts both legs at nonzero hard momenta in this domain.
            assert not np.any((hh<=hi)&(denom<=0))
            I=np.divide((2*np.pi**2/(1,4)[j])*eligible,denom,
                        out=np.zeros((N,N)),where=allowed)
            legs.append(I)
        total=C*np.sum(legs[0]*legs[1]*(E[0,1:][:,None]+E[1,1:][None,:]))/(H*e.sum())
        tests.append(float(total))
    native=high(s,p_cut,16)[1]/(H*e.sum())
    assert 0<native<tests[1]<=tests[0]
    return dict(T_GeV=s["T_visible_GeV"],physical_parent_cut_GeV=p_cut,
                frozen_reconstructed_inverse_energy_over_rho_H=float(native),
                number_only_upper=tests[0],number_and_energy_only_upper=tests[1])


def run(eos,D):
    a=read(D/"warm-light-occupation-next-q32-tape-2026-11-26.json.gz")[0]["snapshot"]["scale_factor"]
    paths=("warm-q64-history-serial-b-2026-11-26.json.gz","warm-history-lower-q64-16-terminal-2026-11-26.json.gz",
           "warm-history-lower-q64-32-terminal-2026-11-26.json.gz")
    out=[]
    for f in paths:
        u=read(D/f)[-1]; s=u["snapshot"]
        assert u["snapshot_sha256"]==digest(s)
        p=one(s,a/s["scale_factor"],eos)
        out.append(dict(source=f,snapshot_sha256=u["snapshot_sha256"],**p))
    return dict(scope="Shape-independent upper for the stored finite q64 retained distributions "
                      "at these checkpoints with exact number or two moments per shell. "
                      "No old exterior occupation, longer gross re-use, additional heavy sources "
                      "or whole-history claim.",rows=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("eos",type=Path)
    p.add_argument("--tapes",type=Path,default=Path("results"));a=p.parse_args()
    print(json.dumps(run(a.eos,a.tapes),indent=2,allow_nan=False))
