"""Reactivate two-packet late exterior examples without changing cohort birth moments."""

import argparse
import json
from pathlib import Path

import numpy as np

from frw_exterior_cohort_reference import Exterior
from frw_exterior_restart_reference import export_restart, import_restart
from joint_prompt_two_moment_evolution_reference import Cells, kin
from late_radiation_background_reference import MPL, TabulatedRadiation
from three_species_decay_onsager_reference import M,m,mp
from warm_exterior_convex_tail_budget_reference import bounds
from warm_light_occupation_stage_audit_reference import digest, read
from warm_prompt_excess_energy_tail_reference import integrate
from warm_prompt_reentry_shape_audit_reference import get_mom,pad


def reactivate(S,hard,half_width,cap=None):
    q=np.asarray(S["lab_momentum_edges_GeV"])
    assert q[-1]<hard-half_width
    wide=np.unique(np.r_[q,np.geomspace(q[-1],hard-half_width,5),hard+half_width,hard*.025+hard])
    if cap is not None:
        assert cap>wide[-1]
        wide=np.r_[wide,np.geomspace(wide[-1],cap,5)[1:]]
    X=import_restart(S["positive_exterior_restart"])
    a,d,j,s,L,n,E=X.array()
    assert np.all(s==1)
    low=L+.20*(E-L)
    high=np.hypot(d,hard/a)
    assert np.all((low<E)&(E<high))
    share=(E-low)/(high-low)
    packet_n=np.array([(1-share)*n,share*n])
    birth=np.array([low,high])
    np.testing.assert_allclose((packet_n*birth).sum(axis=0),n*E,rtol=4e-15)
    pp=a*np.sqrt((birth-d)*(birth+d))
    ix=np.searchsorted(wide,pp,side="right")-1
    assert np.min(ix)>=0 and np.max(ix)<len(wide)-1
    th=pad(S,wide)
    N0,K0,E0=get_mom(wide,th)
    err=[]
    for k,mass in enumerate((mp,m)):
        sel=d==mass
        ids=ix[:,sel].ravel()
        w=packet_n[:,sel].ravel()
        Nadd=np.bincount(ids,w,minlength=len(wide)-1)
        Kadd=np.bincount(ids,w*kin(pp[:,sel].ravel(),mass),minlength=len(wide)-1)
        take=np.flatnonzero(Nadd>0.)
        cell=Cells(wide,mass,(1,4)[k])
        Q,e,_=cell.reconstruct([N0[k]+Nadd,K0[k]+Kadd])
        th[k][:,take]=np.asarray(Q)[:,take]
        err.append(e)
    N,K,F=get_mom(wide,th)
    total_particles=N0.sum()+n.sum()
    direct=float(np.sum(packet_n*np.hypot(pp,d)))
    np.testing.assert_allclose(N.sum(),total_particles,rtol=5e-14)
    np.testing.assert_allclose(F.sum(),E0.sum()+direct,rtol=5e-14)
    assert X.energy(1.)<=direct*(1+2e-14) and direct<=X.upper_bound(1.)*(1+2e-14)
    out=dict(S,lab_momentum_edges_GeV=wide.tolist(),entropy_variable_scalar=th[0].tolist(),
             entropy_variable_fermion=th[1].tolist(),positive_exterior_restart=export_restart(Exterior()),
             represented_species_count_a3=(N.sum(axis=1)*S["scale_factor"]**3).tolist())
    return out, direct, err, pp, packet_n


def run(eos,prior,archive,half_width):
    x=read(archive)
    S=x["snapshot"]
    assert x["snapshot_sha256"]==digest(S)
    T=TabulatedRadiation(eos,S["T_visible_GeV"])
    cut=np.hypot(M,1./read(prior)["branches"][-1]["entropy_expansion_from_original_anchor"])
    E0=get_mom(np.asarray(S["lab_momentum_edges_GeV"]),
                                         [S["entropy_variable_"+k] for k in ("scalar","fermion")])[2].sum()
    ext=import_restart(S["positive_exterior_restart"])
    G=ext.array()
    ceiling=bounds(G,cut).sum(axis=1)
    results=[]
    for hard in (.68,.80,1.):
        A,Es,err,pp,nn=reactivate(S,hard,half_width)
        E=np.hypot(pp,G[1])
        actual=[np.sum(nn*np.maximum(E-cut,0.)),np.sum(nn*E*(E>=cut))]
        assert np.max(np.asarray(actual)/ceiling[1:])<=1+1e-12
        rho=E0+Es
        H=np.sqrt((T.state(1.)[1]+rho)/3)/MPL
        out=[(integrate(A,cut,o,False)/(rho*H)).tolist() for o in (6,10,14)]
        results.append(dict(hard_momentum_GeV=hard,reactivated_snapshot_sha256=digest(A),
                            source_bounds_excess_and_above_cut_GeV4=ceiling[1:].tolist(),
                            point_packet_excess_and_above_cut_GeV4=actual,
                            external_energy_over_old_retained=Es/E0,per_species_moment_recovery=err,
                            positive_and_negative_prompt_ledger_by_order=out))
    return dict(scope="Three frozen positive packet choices from the known ballistic receipts, not "
                      "a shape-independent reacting bound or a solved late kinetic history.",
                base_snapshot_sha256=x["snapshot_sha256"],cut_GeV=cut,hard_shell_half_width_GeV=half_width,quadrature=[6,10,14],cases=results)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","prior","archive"):
        p.add_argument(k,type=Path)
    p.add_argument("--half-width",type=float,default=.005)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.prior,a.archive,a.half_width),indent=2,allow_nan=False))
