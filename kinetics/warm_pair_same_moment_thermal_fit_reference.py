"""Frozen common pair response under separate same-number/energy thermal fits."""

import argparse
import json
from pathlib import Path

import numpy as np

from frw_exterior_restart_reference import import_restart
from joint_prompt_two_moment_evolution_reference import Cells,kin
from joint_visible_heavy_adiabatic_reference import snapshot_functions
from late_radiation_background_reference import MPL,TabulatedRadiation
from three_species_decay_onsager_reference import m,mp
from warm_history_lower_matched_audit_reference import measures
from warm_light_occupation_stage_audit_reference import digest,read
from warm_pair_common_fermion_chart_reference import rate
from warm_prompt_reentry_shape_audit_reference import get_mom


def fit(eos,path):
    v=read(path)
    v=v[-1] if isinstance(v,list) else v
    S=v["snapshot"]; assert digest(S)==v["snapshot_sha256"]
    p=np.asarray(S["lab_momentum_edges_GeV"])
    keys=["entropy_variable_scalar","entropy_variable_fermion"]
    th=[np.asarray(S[k]) for k in keys]
    N,K,E=get_mom(p,th)
    th2=[]
    report=[]
    for j,mass in enumerate((mp,m)):
        c=Cells(p[[0,-1]],mass,(1,4)[j])
        z,err,_=c.reconstruct(np.array([N[j].sum(),K[j].sum()])[:,None])
        z=np.asarray(z)
        slope=z[1,0]/c.dk[0]
        assert slope>0. and err<5e-12
        k=kin(p,mass)
        th2.append(np.stack((z[0,0]+slope*k[:-1],slope*np.diff(k))))
        report.append(dict(temperature_GeV=float(1/slope),mu_GeV=float(mass-z[0,0]/slope),
                           cell_fit_quad_residual=err))
    F=dict(S,**dict(zip(keys,(x.tolist() for x in th2),strict=True)))
    b=measures(S); c=measures(F)
    np.testing.assert_allclose(b[:2].sum(axis=-1),c[:2].sum(axis=-1),rtol=6e-12,atol=0.)
    discrepancy=(abs(b-c).sum(axis=-1)/b.sum(axis=-1)).tolist()
    T=S["T_visible_GeV"]
    extra=import_restart(S["positive_exterior_restart"]).energy(1.)
    rho=TabulatedRadiation(eos,T).state(1.)[1]
    H=np.sqrt((rho+E.sum()+extra)/3)/MPL
    runs=[]
    for order in ((3,7,6,2),(4,10,10,4)):
        u=[]
        for Q in (S,F):
            B,_,L,_=snapshot_functions(Q)
            _,_,g,closure,_=rate(p,L,B,order)
            d=g[0]-g[1]
            u.append(dict(gross_inverse_over_forward_number_energy=(g[0]/g[1]).tolist(),
                         net_fermion_number_over_its_NH=float(2*d[0]/(N[1].sum()*H)),
                         net_fermion_energy_over_total_rhoH=float(d[1]/(E.sum()*H)),ledger=closure))
        runs.append(dict(order=order,actual=u[0],same_thermal_moments=u[1]))
    return dict(T_GeV=T,snapshot_sha256=digest(S),reference_fits=report,
                shell_L1_number_energy_pressure_relative_to_actual=discrepancy,controls=runs)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("paths",nargs="+",type=Path)
    a=p.parse_args()
    print(json.dumps([fit(a.eos,q) for q in a.paths],indent=2,allow_nan=False))
