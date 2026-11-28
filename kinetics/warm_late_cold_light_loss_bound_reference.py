"""Inverse loss upper for a known hard-light count, without a light shape assumption."""

import argparse
import json
from pathlib import Path

import numpy as np

from decay_scalar_ir_relaxation_reference import limits
from frw_exterior_restart_reference import import_restart
from late_radiation_background_reference import MPL,TabulatedRadiation
from three_species_decay_onsager_reference import M,P,mp,m,width
from warm_late_cold_pointwise_bose_dose_reference import envelope
from warm_late_cold_scalar_single_use_receipt_reference import eligible
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom


def inverse_bound(p, archive, parent_cap, n):
    start=envelope(p,parent_cap,n)["cold_reach_floor_GeV"]
    if start>=p:
        return 0.
    # Within one panel, the recoil threshold decreases with p, while
    # the rate prefactor and number-to-unweighted phase factor increase
    # toward different endpoints. Use independent safe endpoints.
    samples=np.geomspace(start,p,n+1)
    E=np.hypot(mp,samples)
    L=np.array([limits(float(x))[1] for x in samples])-E
    g=4/(2*np.pi**2)
    pref=4*M*M/(2*P*samples*E)
    weights=g*L*np.sqrt((L-m)*(L+m))
    high=np.array([eligible(x,archive)[1] for x in L])
    assert np.all(np.diff(L)<0) and np.all(np.diff(high)>=-1e-40)
    return float(np.max(pref[:-1]*high[1:]/weights[1:]))


def run(eos,paths):
    rows=[]
    for p in paths:
        x=read(p)
        S=x["snapshot"]
        assert x["snapshot_sha256"]==digest(S)
        T=S["T_visible_GeV"]
        cut=.1*T
        ep=np.asarray(S["lab_momentum_edges_GeV"])
        _,_,E=get_mom(ep,[S["entropy_variable_"+n] for n in ("scalar","fermion")])
        H=np.sqrt((TabulatedRadiation(eos,T).state(1.)[1]+E.sum())/3.)/MPL
        z=import_restart(S["positive_exterior_restart"]).array()
        thresh=limits(cut)[1]-np.hypot(mp,cut)
        assert thresh>np.hypot(m,ep[-1])
        studies=[]
        for cap in (1.05,2.,3.5,7.5,20.):
            Ecap=float(np.hypot(cap,mp)+np.hypot(cap,m))
            W=[inverse_bound(cut,z,Ecap,n)*width/H for n in (24,96,384)]
            assert W[0]>=W[1]>=W[2]>=0.
            studies.append(dict(common_daughter_cap_GeV=cap,parent_energy_cap_GeV=Ecap,
                                inverse_loss_log_derivative_bound_over_H=W,
                                loss_fraction_upper_for_Hdt_1e_minus7=float(-np.expm1(-W[-1]*1e-7))))
        rows.append(dict(archive=x["snapshot_sha256"],untracked_light_excluded=True,cut_p_over_T=.1,
                         minimum_required_light_energy_GeV=float(thresh),studies=studies))
    return dict(scope="Instantaneous conditional light-loss bound at original stored light moments, with no untracked "
                      "input and finite parent cap. A finite time estimate requires sustaining these tail "
                      "bounds and cap, and pertains only to this inverse channel; no time persistence inferred.",
                cases=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("archives",nargs="+",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.archives),indent=2,allow_nan=False))
