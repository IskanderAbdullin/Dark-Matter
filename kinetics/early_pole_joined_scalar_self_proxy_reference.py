"""Reciprocal scalar decay and a conservative held scalar-scattering proxy."""
import argparse
import json

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import roots_legendre

from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_live_aperture_witness_reprojection_reference import initial_projection
from early_pole_native_aperture_witness_reference import PROFILES
from early_pole_scalar_cold_screened_flux_reference import Screen
from early_pole_scalar_elastic_stencil_reference import Elastic
from early_pole_scalar_hartree_discrete_gradient_reference import gradient


def run(grid,bins,order,duration,steps,subdivision,mixes=(0.,.97),refresh=False):
    if bins<4 or steps<1 or order<2 or not duration>0.:
        raise ValueError("Invalid scalar mesh, time or order.")
    if not mixes or any(not 0<=m<1 for m in mixes):
        raise ValueError("Witness mixtures must lie in [0,1).")
    profile=next(x for x in json.loads(PROFILES.read_text()) if x['grid']==grid)
    C,archive,Y0,Y1,*_=initial_projection(profile,subdivision,bins,duration)
    block=slice(C.j0,C.j0+C.n)
    N0=C.ns.copy()
    _,_,_,f0,m0=C.scalar(N0)
    edges=np.geomspace(1e-6,1e-5,bins+1)
    H=np.sqrt(np.pi**2*106.75/90)*(86./archive['x_end'])/2.435e18
    Q=Elastic(C.q,m0,edges,order,C.lam,Screen(edges,f0,m0,C.lam,order))
    np.testing.assert_array_equal(Q.vol,C.sw)
    tau=duration/steps
    def decay(_,z):
        return Cold.rhs(C,z,False)

    def total(z):
        Y=z[:C.j0].reshape(-1,2)
        _,_,_,e,s=C.M.weights(Y,True)
        _,w,b,_,mass=C.scalar(z[block])
        counts=[float(Y[rr,0].sum()) for rr in C.M.r]
        return e+w,s+b,counts,mass

    calls=0
    self_count_defect=0.
    max_local_self_work=0.
    max_secant_defect=0.
    q_all=np.r_[C.q,C.spq]
    def self_half(z,dt):
        nonlocal calls,self_count_defect,max_local_self_work,max_secant_defect
        start=z[block].copy()
        kernel=Q
        if refresh:
            _,_,_,f,m=C.scalar(start)
            kernel=Elastic(C.q,m,edges,order,C.lam,Screen(edges,f,m,C.lam,order))
        def update(_,logn):
            N=start*np.exp(logn)
            r,_=kernel.flux(N/kernel.vol)
            return r/(H*N)
        q0=np.zeros(len(start))
        span=min(dt/100,.01/max(1.,np.max(abs(update(0.,q0)))))
        h=solve_ivp(update,(0.,dt),q0,method='Radau',
                    rtol=2e-10,atol=2e-12,max_step=dt/8,first_step=span,dense_output=True)
        if not h.success:
            raise RuntimeError(h.message)
        def ledger(order):
            x,w=roots_legendre(order)
            d=(h.t[1:]-h.t[:-1])[:,None]/2
            t=(h.t[1:]+h.t[:-1])[:,None]/2+d*x
            N=start[:,None]*np.exp(h.sol(t.ravel()))
            out=np.zeros(2)
            for n,weight in zip(N.T,(d*w).ravel(),strict=True):
                r,entropy=kernel.flux(n/kernel.vol)
                E,*_=C.scalar(n)
                # Exact Hartree-gradient work on the independently
                # diagnosed fixed-mass scattering graph.
                work=(E-E[0])@r+E[0]*np.sum(r)
                out+=weight*np.r_[entropy,work]/H
            return out
        lo,hi=ledger(3),ledger(5)
        assert np.max(abs(lo-hi))<1e-17
        end=start*np.exp(h.y[:,-1])
        _,r0,_,_,m0=C.scalar(start)
        _,r1,_,_,m1=C.scalar(end)
        secant=gradient(q_all,np.r_[start,C.spect],np.r_[end,C.spect],m0,m1,C.lam)
        exact_work=float(secant[:C.n]@(end-start))
        dead=abs((r1-r0)-exact_work)
        max_secant_defect=max(max_secant_defect,dead)
        assert dead<5e-19 and np.all(np.diff(secant[:C.n])>0)
        self_count_defect=max(self_count_defect,abs(float(end.sum()-start.sum())))
        max_local_self_work=max(max_local_self_work,abs(hi[-1]))
        z[block]=end
        z[-2:]+=hi
        calls+=h.nfev
        return z

    cases=[]
    for mix in mixes:
        for use in (False,True):
            calls=0
            self_count_defect=0.
            max_local_self_work=0.
            max_secant_defect=0.
            z0=np.r_[((1-mix)*Y0+mix*Y1).ravel(),N0,0.,0.]
            r0,s0,c0,_=total(z0)
            z=z0.copy()
            trajectory=[]
            def advance(z,dt):
                D=solve_ivp(decay,(0.,dt),z,method='DOP853',rtol=2e-10,atol=3e-13)
                if not D.success:
                    raise RuntimeError(D.message)
                return D.y[:,-1].copy()
            for k in range(steps+1):
                if k:
                    if use:
                        z=self_half(z,tau/2)
                    z=advance(z,tau)
                    if use:
                        z=self_half(z,tau/2)
                R,S,c,m=total(z)
                gain=float((z[block]-N0).sum())
                cc=np.asarray(c)-c0
                count_defects=np.r_[cc[0],cc[1]+gain,cc[2]-gain]
                assert np.max(abs(count_defects))<5e-12
                assert abs(R-r0)<3e-12 and abs(R-r0-z[-1])<3e-12 and abs(S-s0-z[-2])<3e-11
                trajectory.append(dict(exposure=k*tau,scalar_gain=gain,energy_change=R-r0,
                                       independent_energy_work=z[-1],entropy_closure=S-s0-z[-2],
                                       independent_entropy_ledger=z[-2],
                                       count_defects=count_defects.tolist(),scalar_pole=m,
                                       scalar_l1=float(np.sum(abs(z[block]-N0))/N0.sum()),
                                       instantaneous_cold_source=float(Cold.rhs(C,z,False)[block].sum())))
            assert min(np.diff([s['independent_entropy_ledger'] for s in trajectory]))>=-1e-18
            cases.append(dict(witness_fraction=mix,scalar_elastic_proxy=use,rows=trajectory,
                              scalar_shape_solver_calls=calls,max_isolated_self_count_error=self_count_defect,
                              max_self_halfstep_dispersion_work=float(max_local_self_work),
                              max_hartree_finite_secant_defect=float(max_secant_defect)))
    return dict(grid=grid,bins=bins,order=order,steps=steps,duration=duration,subdivision=subdivision,
                shape_closed_domain_fraction=1-Q.raw_escape,cases=cases,
                **(dict(refreshed_scalar_vertex_and_mass=True) if refresh else {}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=('fine','third'))
    p.add_argument('--bins',type=int,default=16)
    p.add_argument('--order',type=int,default=6)
    p.add_argument('--duration',type=float,default=2e-7)
    p.add_argument('--steps',type=int,default=8)
    p.add_argument('--subdivision',type=int,default=3)
    p.add_argument('--mixes',type=float,nargs='+',default=[0.,.97])
    p.add_argument('--refresh',action='store_true')
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
