"""Four-level time ladder for the coupled q32 support benchmark."""

import argparse
import json
from pathlib import Path

import numpy as np

from frw_exterior_restart_reference import import_restart
from joint_prompt_two_moment_evolution_reference import Cells
from scipy.special import expit
from three_species_decay_onsager_reference import m,mp
from warm_coupled_serial_restart_reference import main,sha
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom


def get(S):
    q=np.asarray(S["lab_momentum_edges_GeV"])
    d=[np.asarray(S["entropy_variable_"+s]) for s in ("scalar","fermion")]
    a=S["scale_factor"]
    N,K,E=get_mom(q,d)
    pressure=[]
    for j,mass in enumerate((mp,m)):
        c=Cells(q,mass,(1,4)[j])
        z=d[j][0,:,None]+d[j][1,:,None]*c.u
        f=expit(-z) if j else 1/np.expm1(z)
        e=mass+c.low[:,None]+c.dk[:,None]*c.u
        p2=(e-mass)*(e+mass)
        pressure.append((c.w*f*p2/(3*e)).sum(axis=1)*c.vol)
    return np.asarray([N,E,pressure])*(a**3)


def audit(eos,source,paths):
    S=read(source)["snapshots"][0]
    q0=np.asarray(S["lab_momentum_edges_GeV"])*S["scale_factor"]
    s=get(S)
    total=s[0].sum(axis=1)
    scale=np.vstack((total,np.tile(s[1].sum(),(2,2))))
    a=import_restart(S["positive_exterior_restart"]).array()
    n0=a.shape[1]
    states=[]; rows=[]
    for n,path in zip((2,4,8,16),paths,strict=True):
        R=read(path); T=R["snapshot"]
        assert R["blocks_finished"]==R["blocks_total"]==n
        assert R["snapshot_sha256"]==digest(T) and R["plan_sha256"]==digest(R["plan"])
        assert len(R["ledger"])==n
        assert R["plan"]["source_sha256"]==sha(source) and R["plan"]["eos_sha256"]==sha(eos)
        q=np.asarray(T["lab_momentum_edges_GeV"])*T["scale_factor"]
        np.testing.assert_allclose(q,q0,rtol=5e-14,atol=0.)
        ext=import_restart(T["positive_exterior_restart"]).array()
        b=T["scale_factor"]/S["scale_factor"]
        assert abs(b-1.001)<5e-14
        c=ext[:,:n0].copy()
        c[[0,5]]*=np.array([b,b**3])[:,None]
        np.testing.assert_allclose(c,a,rtol=5e-14,atol=0.)
        ex=ext[:,n0:]
        exterior_number=np.array([ex[5,ex[1]==m0].sum() for m0 in (mp,m)])*T["scale_factor"]**3
        new_birth=(ex[5]*ex[6]).sum()*T["scale_factor"]**3
        mm=get(T)
        assert abs((mm[0].sum()+exterior_number.sum())/total.sum()-1)<2e-11
        peak0=float(expit(-T["entropy_variable_fermion"][0][0]))
        states.append((mm/scale[:,:,None],new_birth))
        rows.append(dict(blocks=n,visible_T_GeV=T["T_visible_GeV"],
                         checkpoint_sha256=R["snapshot_sha256"],
                         latest_light_rest_occupation=peak0,
                         total_represented_count_change_per_initial_species=(mm[0].sum(axis=1)/total-1).tolist(),
                         exact_new_collisionless_birth_energy_over_initial_dark_total=float(new_birth/s[1].sum()),
                         new_exterior_species_count_per_initial=(exterior_number/total).tolist()))
    errors=[np.sum(abs(states[i+1][0]-states[i][0]),axis=2) for i in range(3)]
    birth=[abs(states[i+1][1]-states[i][1])/max(states[i+1][1],1e-280) for i in range(3)]
    assert np.all(np.asarray(errors)>0.)
    return dict(scope="Matched selected-operator q32 continuation to the same short expansion factor. "
                      "A test of time-grid sensitivity, not the continuum temporal or momentum limit.",
                order_of_L1_metrics=["number","full_energy","pressure"],cases=rows,
                successive_L1_errors=[x.tolist() for x in errors],
                ratios_of_successive_error_levels=[(l/u).tolist() for l,u in zip(errors[:-1],errors[1:],strict=True)],
                successive_relative_change_in_new_exterior_birth_energy=birth)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("source",type=Path)
    p.add_argument("paths",type=Path,nargs=4)
    p.add_argument("--execute",action="store_true")
    a=p.parse_args()
    if a.execute:
        for n,f in zip((2,4,8,16),a.paths,strict=True):
            assert f.name=="checkpoint.json"
            main(a.eos,a.source,f.parent,1.001,n,True,None)
    print(json.dumps(audit(a.eos,a.source,a.paths),indent=2,allow_nan=False))
