"""A reciprocal, live warm spectator extension of the conditional decay-only model."""
import argparse
import gzip
import json

import numpy as np
from scipy.integrate import solve_ivp

from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_joined_scalar_exterior_source_reference import band_kernel
from early_pole_native_aperture_witness_reference import GRID


def run(path,projection,duration,panels,radial_order,outer,linear_probe=False):
    with gzip.open(path,'rt') as f:
        raw=json.load(f)['branches'][projection]
    grid=raw['grid']
    with gzip.open(GRID/f'early-pole-spectrum-{grid}-m06-2026-11-26.json.gz','rt') as f:
        archive=json.load(f)
    shot=raw['snapshots'][-1]
    Y0=np.asarray(shot['hard_moments'])
    N0=np.asarray(shot['cold_bins'])
    n=len(N0)
    C=Cold(grid,1,0,n,True,exposure=raw['duration'],reference_aperture_grid=raw['aperture_grid'],
           archive_projection=projection)
    for x,y in zip(C.M.edge_sets,raw['hard_edge_sets'],strict=True):
        np.testing.assert_array_equal(x,y)
    cold=(C.q.copy(),C.sw.copy())
    spectator=C.spect.copy()
    spine=C.spq
    native_edges=np.asarray(archive['scalar_momentum_edges'])
    volumes=np.diff(native_edges**3)/(6*np.pi**2)
    lo,hi=outer
    take=(spine>=lo)&(spine<=hi)&(spectator>0)
    assert take.any()
    idx=np.flatnonzero(take)
    overlap=np.maximum(0.,np.minimum(native_edges[1:]**3,1e-15)-
                          np.maximum(native_edges[:-1]**3,1e-18))/(6*np.pi**2)
    # The part already represented by live cold bins must not enter the
    # promoted volume twice if a native shell crosses the cold boundary.
    warm_q,warm_w,warm_N=spine[take],(volumes-overlap)[take],spectator[take]
    assert np.all(warm_w>0.)
    initial_pole=C.scalar(N0)[-1]
    initial_energy=C.scalar(N0)[1]

    def evolution(tail):
        monitor=[]
        if tail:
            C.q=np.r_[cold[0],warm_q]
            C.sw=np.r_[cold[1],warm_w]
            C.spect=np.where(take,0.,spectator)
            start=np.r_[N0,warm_N]
        else:
            C.q,C.sw=cold
            C.spect=spectator.copy()
            start=N0.copy()
        C.n=len(start)
        base=np.r_[Y0.ravel(),start,0.,0.]
        _,rho,_,_,m0=C.scalar(start)
        count_defect=float(start.sum()+C.spect.sum()-(N0.sum()+spectator.sum()))
        assert abs(count_defect)<1e-16 and abs(rho-initial_energy)<1e-19 and abs(m0-initial_pole)<4e-18

        def totals(flat):
            Y=flat[:C.j0].reshape(-1,2)
            _,_,_,e,s=C.M.weights(Y,True)
            _,sc,ss,_,p=C.scalar(flat[C.j0:-2])
            return e+sc,s+ss,p,np.array([Y[z,0].sum() for z in C.M.r[1:]])

        def evaluate(_,flat):
            Y=flat[:C.j0].reshape(-1,2)
            N=flat[C.j0:-2]
            Ec,F,Eg,*_=C.M.weights(Y,True)
            es,_,_,fa,pole=C.scalar(N)
            flux=np.zeros_like(Y)
            dN=np.zeros_like(N)
            entropy=0.
            hard_bands=np.zeros((2,2))
            for i,(q,w,f) in enumerate(zip(C.q,C.sw,fa,strict=True)):
                v,d,s=band_kernel(C,Ec,F,q,pole,f,radial_order=radial_order,ledger=True)
                flux+=w*d
                dN[i]=w*v[2]
                entropy+=w*s
                for j,u in enumerate(C.M.r[1:]):
                    hard_bands[int(i>=n),j]+=w*np.sum(Eg[u]*d[u])
            work=np.sum(Eg*flux)+es@dN
            assert abs(work)<5e-13 and entropy>=-1e-30
            monitor.append(float(abs(work)))
            return np.r_[flux.ravel(),dN,entropy,work],hard_bands

        init,_=evaluate(0.,base)
        if not tail:
            np.testing.assert_allclose(init[C.j0:C.j0+n],shot['instantaneous_source_by_q'],
                                       rtol=5e-9,atol=2e-12)
        atol=np.r_[np.full(C.j0,3e-15),np.full(n,2e-18),
                   np.full(len(start)-n,2e-23),3e-16,3e-20]
        sample=np.linspace(0.,duration,9)
        sol=solve_ivp(lambda t,x:evaluate(t,x)[0],(0.,duration),base,rtol=2e-11,atol=atol,
                      max_step=duration/panels,t_eval=sample,method='DOP853')
        assert sol.success
        e0,s0,p0,c0=totals(base)
        e1,s1,p1,c1=totals(sol.y[:,-1])
        change=sol.y[C.j0:-2,-1]-start
        growth=float(change.sum())
        charge=[float(c1[0]-c0[0]+growth),float(c1[1]-c0[1]-growth)]
        assert abs(e1-e0)<2e-13 and abs(e1-e0-sol.y[-1,-1])<2e-13
        assert abs(s1-s0-sol.y[-2,-1])<3e-13 and max(map(abs,charge))<2e-13
        final_flux,hard_bands=evaluate(duration,sol.y[:,-1])
        return dict(warm_active=tail,initial_number_partition_defect=count_defect,initial_gap=p0,final_gap=p1,
                    initial_count_rates=[float(init[C.j0:C.j0+n].sum()),float(init[C.j0+n:-2].sum())],
                    final_count_rates=[float(final_flux[C.j0:C.j0+n].sum()),float(final_flux[C.j0+n:-2].sum())],
                    cold_gain=float(change[:n].sum()),warm_gain=float(change[n:].sum()),
                    warm_count_gains_by_q=change[n:].tolist(),charge_count_defects=charge,
                    energy_closure=e1-e0,entropy_closure=s1-s0-sol.y[-2,-1],
                    maximum_per_evaluation_energy_defect=max(monitor),
                    cold_gain_history=(sol.y[C.j0:C.j0+n].sum(axis=0)-N0.sum()).tolist(),
                    warm_energy_rates_to_parent_and_daughter=hard_bands[1].tolist())

    no=evolution(False)
    yes=evolution(True)
    np.testing.assert_allclose(no['initial_count_rates'][0],yes['initial_count_rates'][0],rtol=2e-13)

    linear={}
    if linear_probe:
        # After the warm evolution arm C is configured for the exact
        # promoted partition. Hold its start state and treat only the
        # additional warm derivative as an infinitesimal direction.
        base_N=np.r_[N0,warm_N]
        _,_,_,fa,pole=C.scalar(base_N)
        E,F,*_=C.M.weights(Y0,True)
        Hflux=np.zeros_like(Y0)
        sflux=np.zeros_like(base_N)
        for j in range(n,len(C.q)):
            a,b,_=band_kernel(C,E,F,C.q[j],pole,fa[j],radial_order=radial_order,ledger=True)
            sflux[j]=C.sw[j]*a[2]
            Hflux+=C.sw[j]*b

        def cold(Y,N):
            E,F,*_=C.M.weights(Y,True)
            _,_,_,f,p=C.scalar(N)
            return sum(w*band_kernel(C,E,F,q,p,z,radial_order=radial_order)[2]
                       for q,z,w in zip(C.q[:n],f[:n],C.sw[:n],strict=True))
        rows=[]
        for offset in (.003,.001,.0003):
            def derivative(a,b):
                return (cold(Y0+a*offset*Hflux,base_N+b*offset*sflux)
                       -cold(Y0-a*offset*Hflux,base_N-b*offset*sflux))/(2*offset)
            u,v,w=(derivative(*b) for b in ((1.,0.),(0.,1.),(1.,1.)))
            rows.append(dict(diagnostic_offset=offset,hard_only=u,gap_only=v,combined=w,
                             linear_additivity_defect=w-u-v))
        linear['initial_added_warm_direction']=rows
    return dict(grid=grid,projection=projection,history_exposure=shot['exposure'],duration=duration,
                time_panels=panels,radial_order=radial_order,warm_node_cuts=[lo,hi],warm_nodes=warm_q.tolist(),
                warm_archive_intervals=np.c_[native_edges[idx],native_edges[idx+1]].tolist(),
                initial_promoted_counts=warm_N.tolist(),fixed_lambda=C.lam,cases=[no,yes],
                **(dict(partial_shell_promoted_volumes=warm_w.tolist()) if np.any(overlap[take]>0.) else {}),
                **linear,scope='Decay and inverse decay only; finite marker support, no optical expansion, other live channels or final abundance.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('path')
    p.add_argument('--projection',choices=('steps','cell-conservative'),default='steps')
    p.add_argument('--duration',type=float,default=2e-6)
    p.add_argument('--panels',type=int,default=4)
    p.add_argument('--radial-order',type=int,default=8)
    p.add_argument('--outer',type=float,nargs=2,default=[.1,5.])
    p.add_argument('--linear-probe',action='store_true')
    a=p.parse_args()
    if a.duration<=0. or a.panels<1 or a.radial_order<2 or not 1e-5<a.outer[0]<a.outer[1]<=5.:
        p.error('Invalid duration, resolution or warm-node limits.')
    print(json.dumps(run(**vars(a)),indent=2,allow_nan=False))
