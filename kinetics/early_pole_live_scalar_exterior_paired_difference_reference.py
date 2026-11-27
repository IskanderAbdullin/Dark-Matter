"""Evolve baseline and the small difference under a shared time grid."""
import argparse
import gzip
import json

import numpy as np

from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_joined_scalar_exterior_source_reference import band_kernel
from early_pole_native_aperture_witness_reference import GRID


def run(path,projection,duration,steps,order,outer,joined=False):
    with gzip.open(path,'rt') as f:raw=json.load(f)['branches'][projection]
    g=raw['grid']
    with gzip.open(GRID/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz','rt') as f:d=json.load(f)
    snap=raw['snapshots'][-1]
    y0=np.asarray(snap['hard_moments'])
    n0=np.asarray(snap['cold_bins'])
    if joined:
        from early_pole_archived_cold_screened_elastic_reference import Screened
        from early_pole_archived_live_compton_reference import ComptonKernel
        from early_pole_archived_live_flavor_conversion_reference import ConversionKernel
        from early_pole_screened_active_response_update_reference import ActiveResponse
        from early_pole_screened_elastic_channel_sum_reference import ElasticChannels
        from early_pole_screened_elastic_discrete_gradient_reference import step
        C=Screened(g,1,len(n0),raw['duration'],(1,1,2),.08,1.,raw['aperture_grid'],projection)
        active=ActiveResponse(C,.12)
        soft=ElasticChannels(C,g,1.)
        add=(ComptonKernel(C.M),ConversionKernel(C.M))
    else:
        C=Cold(g,1,0,len(n0),True,exposure=raw['duration'],reference_aperture_grid=raw['aperture_grid'],
               archive_projection=projection)
    for a,b in zip(C.M.edge_sets,raw['hard_edge_sets'],strict=True):
        np.testing.assert_array_equal(a,b)
    c_q,c_w=C.q.copy(),C.sw.copy()
    orig=C.spect.copy()
    warm=(C.spq>=outer[0])&(C.spq<=outer[1])&(orig>0)
    q=C.spq[warm]; edges=np.array(d['scalar_momentum_edges'])
    overlap=np.maximum(0.,np.minimum(edges[1:]**3,1e-15)-np.maximum(edges[:-1]**3,1e-18))
    vol=((np.diff(edges**3)-overlap)/(6*np.pi**2))[warm]
    assert len(vol)>0 and min(vol)>0.
    w0=orig[warm]
    specs=(orig,np.where(warm,0.,orig))
    banks=((c_q,c_w),(np.r_[c_q,q],np.r_[c_w,vol]))
    h0,nw=len(n0)+C.j0,len(q)
    initial=np.r_[y0.ravel(),n0]
    U=np.zeros(2*h0+nw)

    def set_bank(arm):
        C.spect=specs[arm]
        C.q,C.sw=banks[arm]
        C.n=len(C.q)

    def split(z,arm):
        set_bank(arm)
        y=z[:C.j0].reshape(-1,2)
        N=z[C.j0:]
        E,F,G,*_=C.M.weights(y,True)
        En,_,_,f,m=C.scalar(N)
        rate=np.zeros((2,h0))
        scalar=np.zeros(C.n)
        for i,w in enumerate(C.sw):
            v,H,_=band_kernel(C,E,F,C.q[i],m,f[i],radial_order=order,ledger=True)
            rate[int(i>=len(n0)),:C.j0]+=w*H.ravel()
            scalar[i]=w*v[2]
        rate[0,C.j0:]=scalar[:len(n0)]
        defect=En@scalar+np.sum(G*rate[:,:C.j0].sum(axis=0).reshape(y.shape))
        assert abs(defect)<3e-13
        return rate,scalar[len(n0):]

    def other_slow(Y):
        out=[C.M.pair_rhs(Y,True)[:-2]]
        state=C.M.weights(Y,True)[:3]
        for source in add:
            d,*_=source.collision(*state)
            out.append(d.ravel())
        return np.asarray(out)

    def evaluate(U):
        base=initial+U[:h0]
        diff=U[h0:2*h0]
        first,_=split(base,0)
        second,w=split(np.r_[base+diff,w0+U[2*h0:]],1)
        delta=(second[0]-first[0])+second[1]
        if joined:
            a0=other_slow(base[:C.j0].reshape(-1,2))
            a1=other_slow((base[:C.j0]+diff[:C.j0]).reshape(-1,2))
            first[0,:C.j0]+=a0.sum(axis=0)
            delta[:C.j0]+=(a1-a0).sum(axis=0)
        return np.r_[first[0],delta,w]

    def total(state,arm):
        set_bank(arm)
        y=state[:C.j0].reshape(-1,2)
        _,_,_,en,s=C.M.weights(y,True)
        _,r,b,*_=C.scalar(state[C.j0:])
        return np.array([en+r,s+b])

    e0=total(initial,0)
    e1=total(np.r_[initial,w0],1)
    assert abs(e0[0]-e1[0])<3e-16
    h=duration/steps
    rows=[]
    elastic_drifts=[[],[]]
    difference_roundoff=[]
    def rk(v,dt):
        A=evaluate(v)
        B=evaluate(v+dt*A/2)
        Cc=evaluate(v+dt*B/2)
        D=evaluate(v+dt*Cc)
        return dt*(A+2*B+2*Cc+D)/6
    for i in range(steps):
        if joined:
            U+=rk(U,h/2)
            hb=initial[:C.j0]+U[:C.j0]
            hd=U[h0:h0+C.j0].copy()
            ha=hb+hd
            moved=[]
            for arm,v in enumerate((hb,ha)):
                yy=v.reshape(-1,2)
                active.refresh(yy,'full')
                end,diag=step(soft,yy,h)
                elastic_drifts[arm].append(diag['rho_drift'])
                moved.append(end.ravel()-v)
            U[:C.j0]+=moved[0]
            U[h0:h0+C.j0]+=moved[1]-moved[0]
            hb_new=initial[:C.j0]+U[:C.j0]
            residual=(hb_new+U[h0:h0+C.j0])-(ha+moved[1])
            difference_roundoff.append(float(np.max(abs(residual))))
            U+=rk(U,h/2)
        else:
            U+=rk(U,h)
        if (i+1)%max(1,steps//4)==0 or i+1==steps:
            rows.append(dict(exposure=h*(i+1),cold_gain=float(U[C.j0:h0].sum()),
                             outer_feedback_on_cold=float(U[h0+C.j0:2*h0].sum()),
                             warm_gain=float(U[2*h0:].sum())))
    base=initial+U[:h0]
    t0=total(base,0)-e0
    t1=total(np.r_[base+U[h0:2*h0],w0+U[2*h0:]],1)-e1
    extra={}
    if joined:
        hard_count=np.sum(y0[:,0])
        count=[float(np.sum(x.reshape(-1,2)[:,0])-hard_count)
               for x in (base[:C.j0],base[:C.j0]+U[h0:h0+C.j0])]
        assert max(abs(t0[0]),abs(t1[0]),*map(abs,count))<3e-12
        assert min(t0[1],t1[1])>=0. and max(difference_roundoff)<2e-16
        extra=dict(joined_macro_steps=steps,slow_rk4_half_steps=2*steps,total_hard_count_defect=count,
                   summed_elastic_energy_defects=[float(sum(x)) for x in elastic_drifts],
                   maximum_elastic_hard_pair_reconstruction_residual=max(difference_roundoff))
    else:
        assert max(abs(t0[0]),abs(t1[0]))<2e-13 and min(t0[1],t1[1])>=0.
    return dict(grid=g,projection=projection,duration=duration,common_rk4_steps=steps,radial_order=order,
                energy_closure=[t0[0],t1[0]],entropy_change=[t0[1],t1[1]],rows=rows,
                **(dict(warm_node_cuts=outer) if outer!=[.1,5.] else {}),
                **extra,scope=('Short fixed-scale all-channel matched split on one grid; explicit difference variable.'
                               if joined else 'Same cold hard box and decay-only physics; explicit difference variable on a shared RK4 grid.'))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('path')
    p.add_argument('--projection',choices=('steps','cell-conservative'),default='steps')
    p.add_argument('--duration',type=float,default=7.2e-5)
    p.add_argument('--steps',type=int,default=24)
    p.add_argument('--order',type=int,default=10)
    p.add_argument('--outer',type=float,nargs=2,default=[.1,5.])
    p.add_argument('--joined',action='store_true')
    a=p.parse_args()
    if a.duration<=0 or a.steps<1 or a.order<2 or not 1e-5<a.outer[0]<a.outer[1]<=5.:
        p.error('Invalid positive exposure, steps or order.')
    print(json.dumps(run(**vars(a)),indent=2,allow_nan=False))
