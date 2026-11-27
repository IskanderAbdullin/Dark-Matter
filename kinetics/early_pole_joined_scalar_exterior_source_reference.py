"""The omitted frozen scalar source on the same live hard snapshots; no tail evolution."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import expit, roots_legendre

from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_decay_mode_clock_reference import G, KAPPA
from early_pole_native_aperture_witness_reference import GRID, occup


def band_kernel(C, E, F, q, ms, fs, radial_order=8, angles=16, ledger=False, *, scale_log=0.):
    P,D = C.M.r[1:]
    ep,ed = E[P],E[D]
    a0,a1=ep.T
    b0,b1=ed.T
    assert min(a1)>0 and min(b1)>0
    es = np.hypot(q,ms)
    left = np.maximum((a0-a1)[:,None], (b0-b1)[None,:]+es)
    right = np.minimum((a0+a1)[:,None], (b0+b1)[None,:]+es)
    iis,jjs = np.where(left<right)
    z,w=roots_legendre(radial_order)
    output=np.zeros(4)
    derivative=np.zeros_like(E) if ledger else None
    entropy=0.
    if ledger and fs<=0.:
        raise ValueError("A differentiable entropy check requires positive scalar occupation.")
    for i,j in zip(iis,jjs,strict=True):
        ia,ib=P.start+i,D.start+j
        A,B=C.M.cells[ia],C.M.cells[ib]
        def ea(p):return a0[i]+a1[i]*A.basis(p)
        def pb(p):return B.pfrombare(B.ec+B.de*(ea(p)-es-b0[j])/b1[j])
        p0,p1=A.pfrombare(A.ec+A.de*(np.array([left[i,j],right[i,j]])-a0[i])/a1[i])
        v=np.linspace(p0,p1,angles+1)
        knots=[p0,p1]
        for cond in (lambda p:pb(p)-p-q,lambda p:pb(p)-p+q,lambda p:pb(p)+p-q):
            y=cond(v)
            knots.extend(v[y==0].tolist())
            for n in np.flatnonzero(y[:-1]*y[1:]<0):
                knots.append(brentq(cond,v[n],v[n+1],xtol=5e-14))
        s=np.unique(knots)
        lo,hi=s[:-1],s[1:]
        pm=(lo+hi)/2
        km=pb(pm)
        take=(abs(km-pm)<=q)&(q<=km+pm)
        if not np.any(take):continue
        h=(hi[take]-lo[take])[:,None]/2
        pp=(pm[take,None]+h*z).ravel()
        dd=pb(pp)
        assert np.all(abs(dd-pp)<=q+2e-12) and np.all(dd+pp>=q-2e-12)
        E0=ea(pp)
        vD=(b1[j]/B.de)*dd/np.hypot(dd,B.m)
        wq=(h*w).ravel()*pp*dd/(E0*(E0-es)*vD)
        ua,ub=A.basis(pp),B.basis(dd)
        a,b=F[ia,0]+F[ia,1]*ua,F[ib,0]+F[ib,1]*ub
        fp,fd=expit(a),expit(b)
        seed=fp*(1-fd)
        net=seed+fs*(fp-fd)
        out=np.array([seed*(1+fs),fd*(1-fp)*fs,net,net*es])
        output+=out@wq
        if ledger:
            n=net*wq
            N=n.sum()
            derivative[ia]-=[N,n@ua]
            derivative[ib]+=[N,n@ub]
            entropy+=n@(a-b+np.log1p(1./fs))
    x=.05 if scale_log==0. else .05*np.exp(scale_log)
    scale=G*KAPPA*x*((1.06*x)**2-ms*ms)/(q*es)
    if ledger:
        return output*scale,derivative*scale,entropy*scale
    return output*scale


def run(path, projection, q_order, radial_order, snapshots):
    with gzip.open(path,'rt') as f:raw=json.load(f)['branches'][projection]
    grid=raw['grid']
    with gzip.open(GRID/f'early-pole-spectrum-{grid}-m06-2026-11-26.json.gz','rt') as f:a=json.load(f)
    C=Cold(grid,1,0,len(raw['snapshots'][0]['cold_bins']),True,exposure=raw['duration'],
           reference_aperture_grid=raw['aperture_grid'],archive_projection=projection)
    assert len(C.M.edge_sets)==len(raw['hard_edge_sets'])
    for x,y in zip(C.M.edge_sets,raw['hard_edge_sets'],strict=True):
        np.testing.assert_array_equal(x,y)
    qref=np.asarray(a['scalar_momentum_edges'])
    counts=np.maximum(0.,np.asarray(a['endpoint_scalar_signed_shell_number_over_T3']))
    fref=counts/(np.diff(qref**3)/(6*np.pi**2))
    cold=np.geomspace(1e-6,1e-5,C.n+1)
    bounds=[1e-8,1e-6,1e-5,.001,.1,5.]
    z,w=roots_legendre(q_order)
    cases=[]
    choices={'first':raw['snapshots'][0],'last':raw['snapshots'][-1]}
    selected = ([(f"step{r['step']}",r) for r in raw['snapshots']] if 'all' in snapshots
                else [(k,choices[k]) for k in snapshots])
    for key,s in selected:
        Y=np.asarray(s['hard_moments'])
        N=np.asarray(s['cold_bins'])
        _,_,_,fs,ms=C.scalar(N)
        E,F,*_=C.M.weights(Y,True)
        central=np.array([band_kernel(C,E,F,q,ms,f,radial_order)[2]
                          for q,f in zip(C.q,fs,strict=True)])*C.sw
        guard=float(np.max(abs(central-np.array(s['instantaneous_source_by_q']))))
        assert guard<1e-10,guard
        eP=E[C.M.r[1]]
        eD=E[C.M.r[2]]
        assert max(eP[:,0]+eP[:,1])-min(eD[:,0]-eD[:,1])<np.hypot(bounds[-1],ms)
        rows=[]
        for lo,hi in zip(bounds[:-1],bounds[1:],strict=True):
            b=np.unique(np.r_[lo,hi,qref[(lo<qref)&(qref<hi)],cold[(lo<cold)&(cold<hi)]])
            # Use log q so both near-IR and broad bands are comparably sampled.
            a0,b0=np.log(b[:-1]), np.log(b[1:])
            h=(b0-a0)[:,None]/2
            pp=np.exp(((a0+b0)[:,None]/2+h*z).ravel())
            ww=(h*w).ravel()*pp**3/(2*np.pi**2)
            f=np.where((cold[0]<=pp)&(pp<=cold[-1]),occup(pp,cold,fs),occup(pp,qref,fref))
            seq=np.array([band_kernel(C,E,F,q,ms,v,radial_order) for q,v in zip(pp,f,strict=True)])
            result=ww@seq
            rows.append(dict(interval=[lo,hi],quad_size=len(pp),gross_gain=float(result[0]),
                             gross_inverse=float(result[1]),net=float(result[2]),
                             scalar_energy_source=float(result[3])))
        cases.append(dict(snapshot=key,step=s['step'],exposure=s['exposure'],scalar_pole=ms,
                          exact_live_cold_node_rate=float(central.sum()),node_replay_defect=guard,rows=rows))
    return dict(grid=grid,projection=projection,q_order=q_order,radial_order=radial_order,cases=cases,
                scope='Frozen tagged scalar source at live snapshots on the same finite hard-leg support; no exterior evolution.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('bundle')
    p.add_argument('--projection',choices=('steps','cell-conservative'),default='steps')
    p.add_argument('--q-order',type=int,default=5)
    p.add_argument('--radial-order',type=int,default=8)
    p.add_argument('--snapshots',nargs='+',choices=('first','last','all'),default=['first','last'])
    a=p.parse_args()
    if min(a.q_order,a.radial_order)<2:
        p.error('Use quadrature orders of at least two.')
    print(json.dumps(run(a.bundle,a.projection,a.q_order,a.radial_order,a.snapshots),
                     indent=2,allow_nan=False))
