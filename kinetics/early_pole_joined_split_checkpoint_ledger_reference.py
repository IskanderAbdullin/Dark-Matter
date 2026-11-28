"""Diagnose a fixed-scale split checkpoint against the unchanged terminal ledgers."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_cold_screened_split_reference import state_digest


def score(d):
    p=d['physics']
    if p.get('witness_mix') or not d['controls'] or d['steps_done']!=p['count']:
        raise ValueError('Use a complete native (not prepared-witness) elastic checkpoint.')
    model=Cold(p['grid'],1,0,p['scalar_bins'],True,exposure=p['duration'],
               reference_aperture_grid=p['aperture_grid'],archive_projection=p['projection'],
               aperture_subdivision=p['aperture_subdivision'])
    z0=np.r_[model.M.init.ravel(),model.ns,0.,0.]
    z=np.asarray(d['state'])
    assert state_digest(z0)==p['state_initial'] and state_digest(z)==d['state_digest']

    def total(z):
        _,_,_,e,s=model.M.weights(z[:model.j0].reshape(-1,2),True)
        _,rho,ent,*_=model.scalar(z[model.j0:model.j0+model.n])
        return e+rho,s+ent

    e0,s0=total(z0)
    e1,s1=total(z)
    err,ent=e1-e0,s1-s0-z[-2]
    drift=[x['rho_drift'] for x in d['controls']]
    running=[math.fsum(drift[:n+1]) for n in range(len(drift))]
    ix=sorted(set([p['count']//4,p['count']//2,2*p['count']//3,3*p['count']//4,p['count']]
                  +[n for n in (300,360) if 0<n<=p['count']]))
    return dict(grid=p['grid'],projection=p['projection'],steps=d['steps_done'],duration=p['duration'],
                final_energy_defect=float(err),final_entropy_defect=float(ent),
                signed_sum_of_elastic_defects=running[-1],energy_gate=3e-12,entropy_gate=3e-11,
                pass_gate=bool(abs(err)<3e-12 and abs(ent)<3e-11),
                largest_per_step_elastic_defect=max(map(abs,drift)),
                first_aggregate_defect_exceeds_energy_gate=next((i+1 for i,v in enumerate(running)
                                                                 if abs(v)>3e-12),None),
                signed_nonelastic_ledger_remainder=float(err-running[-1]),
                provisional_gain=float(d['states'][-1]-d['states'][0]),
                figure_of_merit_at_selected_steps={i:running[i-1] for i in ix if i>0})


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('checkpoint')
    a=p.parse_args()
    print(json.dumps(score(json.loads(Path(a.checkpoint).read_text())),indent=2,allow_nan=False))
