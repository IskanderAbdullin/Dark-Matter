"""Resolve short combined-versus-timelike effects at the exact hot duration."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_tabulated_entropy_split_expansion import Stage
from sm_dark_entropy_exchange_reference import StandardModel


def timelike(eos,p,n,a,dt):
    stage=Stage(eos,p,n)
    S0,E0=stage.entropy(0.,a),stage.energy(0.,a)

    def rhs(t,z):
        x=z[0]
        g,v,ep=stage.rhs(x,z[1:4])
        P=eos.state(z[1])[1]+stage.bath(z[2:4])['pressure']
        return dt*np.r_[g,v,ep/S0,3*g*np.exp(3*x)*P/E0]

    solve=solve_ivp(rhs,(0.,1.),np.r_[0.,a,0.,0.],method='DOP853',rtol=3e-12,atol=3e-14)
    if not solve.success:
        raise ArithmeticError('No exact-duration continuous timelike comparison')
    ref=solve.y[:,-1]
    if max(abs(stage.energy(ref[0],ref[1:4])/E0-1.+ref[5]),
           abs(stage.entropy(ref[0],ref[1:4])/S0-1.-ref[4]))>2e-11:
        raise ArithmeticError('Continuous duration comparison failed either ledger')
    verified=[]
    for k in (4,8):
        x=0.
        u=a.copy()
        w,e=0.,0.
        for _ in range(k):
            x,u,h,_,j,_=stage.step(x,u,dt/k)
            e+=h
            w+=j
        defects=[abs(stage.entropy(x,u)/S0-1.-e/S0),abs(stage.energy(x,u)/E0-1.+w/E0)]
        if max(defects)>1e-11 or min(e,w)<=0.:
            raise ArithmeticError('Adiabatic and collision proof failed')
        verified.append(dict(steps=k,state_and_normalized_budgets=np.r_[x,u,e/S0,w/E0].tolist(),ledger_defects=defects))
    return ref,verified


def run(path,lock):
    raw=path.read_bytes()
    rec=lock.read_bytes()
    d=json.loads(rec)
    digest=hashlib.sha256(raw).hexdigest()
    if d['table_input_sha256']!=digest or d['timelike_comparison']['samples'][-1]['same_duration']:
        raise ArithmeticError('Expected the exact unresolved hot-duration comparison')
    z=d['rows'][-1]
    a=np.asarray(z['initial_state'])
    if a[0]!=.6882487043:
        raise ArithmeticError('Wrong numerical state')
    dt=z['horizon']
    eos=StandardModel(path)
    p,q=[PositivePortal(o) for o in ((13,36,26),(17,44,32))]
    (u,v),(s,t)=[timelike(eos,x,n,a,dt) for x,n in ((p,112),(q,168))]
    a1=np.asarray(z['independent_terminal_coordinates'])
    shifts=a1-s
    if np.max(abs(u-s))>3e-8:
        raise ArithmeticError('Independent timelike quadrature did not reproduce the endpoint')
    diff=np.array([z['levels'][-1]['log_a'],*z['levels'][-1]['state'],
                 z['levels'][-1]['relative_comoving_entropy_added'],z['levels'][-1]['normalized_pressure_work']])-t[-1]['state_and_normalized_budgets']
    return dict(fixed_initial_state=a.tolist(),horizon=dt,table_input_sha256=digest,combined_reference_sha256=hashlib.sha256(rec).hexdigest(),
                fields=['log_a','T_vis','eta_dark','ln_T_dark','comoving_entropy_added_normalized','pressure_work_normalized'],
                combined_continuous=z['independent_terminal_coordinates'],timelike_continuous=s.tolist(),
                continuous_combined_minus_timelike=shifts.tolist(),
                eight_step_combined_minus_timelike=diff.tolist(),timelike_steps=t,
                max_timelike_independent_mesh_difference=np.max(abs(u-s)),
                scope='Direct and split 0.688-GeV scenarios compared at identical physical duration, fixed numerical model and initial state, not cosmological abundance.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('eos',type=Path)
    p.add_argument('combined',type=Path)
    o=p.parse_args()
    print(json.dumps(run(o.eos,o.combined),indent=2,allow_nan=False))
