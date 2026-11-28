"""Reintegrate instantaneous cold source on fixed, already evolved hard snapshots."""
import argparse
import json
from math import fsum
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.special import expit, roots_legendre

from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_decay_mode_clock_reference import G,KAPPA
from early_pole_daughter_direct_cold_volume_quadrature import angular_edge_modes
from early_pole_daughter_local_fraction_curvature_attribution import initial
from early_pole_live_pauli_gauge_moment_decay_reference import Tool
from early_pole_joined_temporal_refinement_reference import conservation


def node_integral(C,Y,n,order,frozen_knots=None):
    M=C.M
    e,f=M.weights(Y,True)[:2]
    _,_,_,_,mass=C.scalar(n)
    x=.05*np.exp(getattr(M,'eta',0.))
    pre=G*KAPPA*x*((1.06*x)**2-mass*mass)/(6*np.pi**2)
    old=np.geomspace(1e-6,1e-5,9)**3
    knots=angular_edge_modes(M,e,mass,257)
    cuts=np.unique(np.r_[old,(knots if frozen_knots is None else frozen_knots)**3])
    z,w=roots_legendre(order)
    rates=[[] for _ in range(8)]
    for a,b in zip(cuts[:-1],cuts[1:],strict=True):
        slot=np.searchsorted(old,(a+b)/2)-1
        B=6*np.pi**2*n[slot]/(old[slot+1]-old[slot])
        h=(b-a)/2
        for vol,wq in zip((a+b)/2+h*z,h*w,strict=True):
            q=np.cbrt(vol)
            W=np.hypot(q,mass)
            model=SimpleNamespace(n=M.decay.n,cells=M.decay.cells,q=q,nod=M.z,w=M.w)
            for i,k,u,v,ww in Tool.events(model,dict(coeff=e[M.start:],es=W)):
                i+=M.start
                k+=M.start
                F=expit(f[i,0]+f[i,1]*u)
                D=expit(f[k,0]+f[k,1]*v)
                rates[slot].append(wq*pre/(q*W)*float(ww@(F*(1-D)+B*(F-D))))
    return dict(gap=float(mass),angular_modes=knots.tolist(),mode_grouped_rate=[fsum(v) for v in rates])


def compare(bundle,native,local,order,hold_startup_cuts=False):
    C,seed,(y0,_,y1),_,depth=initial(bundle)
    if depth!=1:
        raise ValueError('The trajectory in this diagnostic has original hard support.')
    control=('projection','duration','grid','aperture_grid','hard_edge_sets','strength','pair','method',
             'same_soft_strength','compton_strength','conversion_strength','screening_update','response_smear')
    old=None
    output={}
    frozen_knots=None
    if hold_startup_cuts:
        m0=C.scalar(seed)[-1]
        frozen_knots=np.unique(np.concatenate([
            angular_edge_modes(C.M,C.M.weights(y,True)[0],m0,257) for y in (y0,y1)]))
    for a,tag,Y0 in ((native,'native',y0),(local,'local-daughter',y1)):
        if (any(a[k]!=native[k] for k in control) or a['steps']!=64
                or a['screening_update']!='frozen' or a.get('prepared_charged_swap','native')!=tag):
            raise ValueError('Wrong original held first-window trajectory.')
        if a['elastic_energy_budget']!=5e-15:
            raise ValueError('Inverse interpolation needs accepted fixed-inner reference.')
        conservation(a)
        p=a['snapshots']
        np.testing.assert_array_equal(p[0]['hard_moments'],Y0)
        np.testing.assert_array_equal(p[0]['cold_bins'],seed)
        for x,y in zip(a['hard_edge_sets'],C.M.edge_sets,strict=True):
            np.testing.assert_array_equal(x,y)
        t=[s['exposure'] for s in p]
        if old is None:
            old=t
        elif t!=old:
            raise ValueError('Mismatched saved time samples.')
        rows=[]
        for u in p:
            Y=np.asarray(u['hard_moments'])
            n=np.asarray(u['cold_bins'])
            dr=Cold.rhs(C,np.r_[Y.ravel(),n,0.,0.],False)[C.j0:-2]
            np.testing.assert_array_equal(dr,u['instantaneous_source_by_q'])
            m=node_integral(C,Y,n,order,frozen_knots)
            m['original_mode_rates']=dr.tolist()
            rows.append(m)
        output[tag]=rows
    return dict(**(dict(frozen_angular_knots=frozen_knots.tolist()) if hold_startup_cuts else {}),
                scope='Postprocess accepted held split trajectories without evolving a new state; dense scalar sampling assumes constant density within each saved phase-volume bin.',
                quadrature_order=order,exposure=old,arms=output)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('bundle')
    p.add_argument('native')
    p.add_argument('local')
    p.add_argument('--order',type=int,default=8)
    p.add_argument("--hold-startup-cuts",action="store_true")
    a=p.parse_args()
    load=lambda f:json.loads(Path(f).read_text())
    print(json.dumps(compare(a.bundle,load(a.native),load(a.local),a.order,a.hold_startup_cuts),indent=2,allow_nan=False))
