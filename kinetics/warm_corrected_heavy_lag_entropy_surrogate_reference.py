"""Driven heavy lag across a substep-sampled daughter entropy surrogate."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from frw_exterior_restart_reference import import_restart
from joint_visible_heavy_adiabatic_reference import parent_branch, snapshot_functions
from late_radiation_background_reference import MPL, TabulatedRadiation
from three_species_decay_onsager_reference import M,width
from warm_heavy_relaxation_lag_reference import update
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom


def run(eos,parent,end):
    A,B=map(read,(parent,end))
    assert A["snapshot_sha256"]==B["parent_snapshot_sha256"]
    S,T=A["snapshot"],B["snapshot"]
    for q,x in ((A,S),(B,T)):
        assert q["snapshot_sha256"]==digest(x)
    a=T["scale_factor"]/S["scale_factor"]
    span=np.log(a)
    p0=np.asarray(S["lab_momentum_edges_GeV"])
    np.testing.assert_allclose(p0/a,T["lab_momentum_edges_GeV"],rtol=5e-14)
    keys=["entropy_variable_"+z for z in ("scalar","fermion")]
    Q0,Q1=[np.stack([x[k] for k in keys]) for x in (S,T)]
    eps=[import_restart(z["positive_exterior_restart"]).energy(1.) for z in (S,T)]
    R=TabulatedRadiation(eos,S["T_visible_GeV"])
    # Integrate the lag directly, not as the difference of two almost
    # equal level values. Verify the algebraic identity before using it.
    for L in (.01,.3,3.,10.):
        x,h,k=.014,.013,.016
        y=(x-h)*np.exp(-L)+(k-h)*np.expm1(-L)/L
        assert abs(y-(update(x,h,k,L)-k))<5e-18
    grid=np.r_[0.,np.geomspace(1e-4,1,72)]
    results=[]
    for steps,order,radial in ((4,8,8),(8,8,8),(16,8,8),(16,12,8),(16,16,8),(16,16,12)):
        z,w=roots_legendre(radial)
        P=(grid[:-1,None]+np.diff(grid)[:,None]*(1+z)/2).ravel()
        W=(np.diff(grid)[:,None]*w/2).ravel()
        times=np.linspace(0,1,steps+1)
        cells=[]
        for u in times:
            b=np.exp(span*u)
            p=p0/b
            th=(1-u)*Q0+u*Q1
            snap=dict(lab_momentum_edges_GeV=p,**dict(zip(keys,th,strict=True)))
            f,x,g,y=snapshot_functions(snap)
            Mtot=get_mom(p,th)[2].sum()
            H=np.sqrt((R.state(b)[1]+Mtot+(1-u)*eps[0]+u*eps[1])/3)/MPL
            E=np.hypot(M,P/b)
            C=parent_branch(E,f,g,x,y,angle=order)
            cells.append((C[:,0],width*M*C[:,-1]/(E*H),E,Mtot))
        rows=[]
        diff=np.zeros_like(cells[0][0])
        for j,(f,r,E,rho) in enumerate(cells[1:],1):
            L=np.minimum(cells[j-1][1],r)*span/steps
            diff=diff*np.exp(-L)+(f-cells[j-1][0])*np.expm1(-L)/L
            b=np.exp(span*times[j])
            measure=2/np.pi**2*E*P**2*W/(b**3*rho)
            rows.append(dict(fraction=times[j],min_integrated_relaxation=float(min(L)),
                             peak_absolute_occupation_lag=float(max(abs(diff))),
                             energy_absolute_lag_over_retained=float(measure@abs(diff)),
                             frozen_heavy_energy_over_retained=float(measure@f)))
        results.append(dict(substeps=steps,angular_order=order,radial_order=radial,rows=rows))
    return dict(scope="Daughter entropy coefficients linearly interpolated in log expansion on moving "
                      "momentum shells. Then equilibrium and relaxation are sampled on this prescription, "
                      "not enclosed between samples. No active heavy feedback, outside, or other sources.",
                parent_sha256=A["snapshot_sha256"],endpoint_sha256=B["snapshot_sha256"],controls=results)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","parent","endpoint"):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.parent,a.endpoint),indent=2,allow_nan=False))
