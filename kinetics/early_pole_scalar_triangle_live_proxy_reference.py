"""Sensitivity of an isolated scalar shape proxy to triangle-compatible remapping."""
import argparse
import json

import numpy as np
from scipy.integrate import solve_ivp

from early_pole_live_aperture_witness_reprojection_reference import initial_projection
from early_pole_native_aperture_witness_reference import PROFILES
from early_pole_scalar_cold_screened_flux_reference import Screen
from early_pole_scalar_elastic_stencil_reference import Elastic
from early_pole_scalar_triangle_secant_fixture_reference import project


def run(grid,bins,order,duration):
    if bins<4 or order<2 or not duration>0:
        raise ValueError("Invalid finite graph controls.")
    p=next(x for x in json.loads(PROFILES.read_text()) if x['grid']==grid)
    C,archive,Y0,Y1,*_=initial_projection(p,3,bins,duration)
    n0=C.ns.copy()
    e0,r0,s0,f0,m0=C.scalar(n0)
    ed=np.geomspace(1e-6,1e-5,bins+1)
    H=np.sqrt(np.pi**2*106.75/90)*(86./archive['x_end'])/2.435e18
    Q=Elastic(C.q,m0,ed,order,C.lam,Screen(ed,f0,m0,C.lam,order),keep_lab=True)
    times=np.r_[0,np.geomspace(duration/128,duration,11)]
    anchor=int(np.argmax(n0))
    free=np.flatnonzero(np.arange(bins)!=anchor)
    total=float(n0.sum())
    def counts(y):
        n=n0.copy()
        n[free]*=np.exp(y)
        n[anchor]=total-np.sum(n[free])
        if n[anchor]<=0.:
            raise ValueError("A solver trial left the positive fixed-number simplex.")
        return n
    cases=[]
    for panel in (None,0,1):
        def rate(N):
            f=N/Q.vol
            R,A=Q.reaction(f)
            ds=np.log1p(1/f)
            if panel is None:
                F=Q.S.T@R
                h=R@A
                return F/H,h/H,1.,0.
            E,*_=C.scalar(N)
            S,keep,geom=project(Q,E,panel)
            mobility=np.divide(R,A,out=np.zeros_like(R),where=A!=0)
            assert np.min(mobility)>=0
            aff=S@ds
            F=S.T@(mobility*aff)
            h=mobility@(aff*aff)
            return F/H,h/H,float(Q.weight[keep].sum()/Q.weight.sum()),geom
        def fun(_,ln):
            n=counts(ln)
            return (rate(n)[0]/n)[free]
        z0=np.zeros(len(free))
        first=min(duration/100,.01/max(1.,np.max(abs(fun(0.,z0)))))
        h=solve_ivp(fun,(0.,duration),z0,method='Radau',t_eval=times,
                     rtol=2e-9,atol=3e-11,first_step=first,max_step=duration/12)
        if not h.success:raise RuntimeError(h.message)
        rows=[]
        for t,u in zip(times,h.y.T,strict=True):
            n=counts(u)
            E,er,es,*_=C.scalar(n)
            rn,entropy,kept,geom=rate(n)
            cold=[]
            for mix in (0.,.97):
                Y=(1-mix)*Y0+mix*Y1
                out=C.rhs(np.r_[Y.ravel(),n,0.,0.],False)[C.j0:C.j0+C.n]
                cold.append(float(out.sum()))
            assert abs(n.sum()-n0.sum())/n0.sum()<1e-8
            assert abs(er-r0)<1e-13 and entropy>=0
            rows.append(dict(exposure=t,scalar_shape_l1=float(np.sum(abs(n-n0))/n0.sum()),
                             scalar_number_change=float(n.sum()-n0.sum()),scalar_hartree_energy_change=er-r0,
                             scalar_energy_rate_residual=float((E-E[0])@rn+E[0]*np.sum(rn)),
                             entropy_change=es-s0,entropy_rate=entropy,
                             retained_initial_weight=kept,max_momentum_triangle_violation=geom,
                             frozen_charge_step_and_witness_sources=cold))
        assert min(np.diff([a['entropy_change'] for a in rows]))>-1e-27
        cases.append(dict(momentum_panel=panel,samples=rows,rhs_evaluations=h.nfev))
    return dict(grid=grid,bins=bins,order=order,exposure=duration,cases=cases,
                coordinates="positive fixed-count simplex",
                scope="Isolated scalar sensitivity with a held illustrative screened graph; not a physical angle correction.")


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=('fine','third'))
    p.add_argument('--bins',type=int,default=16)
    p.add_argument('--order',type=int,default=6)
    p.add_argument('--duration',type=float,default=2e-7)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
