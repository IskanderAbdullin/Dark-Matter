"""Frozen-mass/frozen-vertex narrow-IR relaxation; a controlled numerical case."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from early_pole_decay_mode_clock_reference import X
from early_pole_massive_daughter_energy_stencil_reference import MassiveFermion
from early_pole_scalar_cold_screened_flux_reference import Screen
from early_pole_scalar_elastic_stencil_reference import Elastic


def run(spec,lo,angular,calculation_scale=1.):
    m=spec['checkpoints'][-1]['scalar_pole_over_T']
    p=np.asarray(spec['scalar_momentum_nodes'])
    edges=np.asarray(spec['scalar_momentum_edges'])
    f=np.maximum(np.asarray(spec['endpoint_scalar_occupations']),0.)
    L=np.searchsorted(edges,lo)
    R=np.searchsorted(edges,1e-5)
    assert edges[R]==1e-5 and R-L>=6
    lm=7.65175e-8
    response=Screen(edges[:R+1],f[:R],m,lm,angular)
    model=Elastic(p[L:R],m,edges[L:R+1],angular,lm,
                  response if calculation_scale==1. else None)
    initial=model.vol*f[L:R]
    count=initial.sum()
    y0=initial/count
    v0,H=model.flux(f[L:R])
    native_initial_clock=np.sum(abs(v0))/count
    assert min(y0)>0.
    def rhs(_,y):
        # Endpoint positivity is checked and reported. Clip only trial
        # states in the numerical Jacobian; tolerated negative counts
        # do not become part of any archived physical solution.
        yp=np.maximum(y,0.)
        v,ent=model.flux(yp*count/model.vol)
        return v/(count*native_initial_clock)
    ts=np.linspace(0,3.,25)
    solution=solve_ivp(rhs,(0.,3.),y0,method='BDF',t_eval=ts,rtol=2e-9,
                       atol=np.maximum(1e-15,y0*1e-11),max_step=.08)
    assert solution.success,solution.message
    Y=solution.y
    assert min(Y.ravel())>=0.
    nr=np.max(abs(Y.sum(axis=0)-1))
    er=np.max(abs((model.K@Y)/(model.K@y0)-1))
    assert max(nr,er)<2e-10
    F=Y*count/model.vol[:,None]
    eta=model.vol@(np.log1p(F)+F*np.log1p(1/F))
    assert min(np.diff(eta))>=-1e-12*eta[-1]
    top=int(np.argmax(f[L:R]))
    h=np.sqrt(np.pi**2*106.75/90)*(86./spec['x_end'])/2.435e18
    frac=count/np.sum(np.diff(edges**3)*f/(6*np.pi**2))
    # Counterfactual source audit: keep parent, fermion daughter and epoch
    # fixed; replace only the isolated scalar spectrum. Recompute the gap.
    control={(67,112,122):(65,90,72),(97,156,176):(95,125,108)}
    dep=MassiveFermion(*control[tuple(spec['cells'])],spec['mass_ratio'],0.,True)
    parent=spec['endpoint_parent_occupations']
    daughter=spec['endpoint_fermion_occupations']
    t=np.log(spec['x_end']/X)
    def coupling(fs):
        y=np.r_[fs,daughter,parent]
        decay=dep.rhs(t,y)[:dep.ns]
        rate=dep.d.w@decay/(dep.d.w@fs)
        gap=dep.d.gap(spec['x_end'],fs)
        return float(rate),float(gap)
    old,gap=coupling(f)
    moved=f.copy()
    moved[L:R]=F[:,-1]
    new,newgap=coupling(moved)
    assert abs(newgap/gap-1)<1e-6
    return dict(x=spec['x_end'],mass_ratio=spec['mass_ratio'],angular_order=angular,
                support_cells=len(y0),support_over_T=edges[[L,R]].tolist(),
                initial_number_fraction_represented=float(frac),
                initial_occupation_peak=float(f[L:R].max()),
                initial_formal_clock_over_Href=float(native_initial_clock/h),
                peak_cell_initial_number_fraction=float(Y[top,0]),
                peak_cell_final_number_fraction=float(Y[top,-1]),
                signed_peak_number_change_relative=float(Y[top,-1]/Y[top,0]-1),
                final_total_variation_number_fraction=float(np.sum(abs(Y[:,-1]-Y[:,0]))/2),
                entropy_fraction_change=float(eta[-1]/eta[0]-1),
                number_closure=float(nr),kinetic_energy_closure=float(er),
                formal_Href_dt_at_end=float(3*h/native_initial_clock),
                dimensionless_evolution_time=3.,screened=calculation_scale==1.,
                initial_decay_net_growth_rate_log_a=old,counterfactual_decay_net_growth_rate_log_a=new,
                updated_gap_ratio=newgap/gap,isolated_support_final_occupation=F[:,-1].tolist())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('spec')
    p.add_argument('--lower',type=float,default=1e-6)
    p.add_argument('--angle',type=int,default=8)
    p.add_argument('--bare',action='store_true')
    a=p.parse_args()
    f=Path(a.spec)
    raw=gzip.decompress(f.read_bytes()) if f.suffix=='.gz' else f.read_bytes()
    d=run(json.loads(raw),a.lower,a.angle,0. if a.bare else 1.)
    print(json.dumps(d,indent=2,allow_nan=False))
