"""Project the moment-matched bare witness into the dressed finite live closure."""
import argparse
import gzip
import json

import numpy as np
from scipy.special import expit, roots_legendre

from early_pole_archived_extended_pair_decay_window_reference import fit_archive_cell
from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_native_aperture_witness_reference import GRID, PROFILES, refine, mode
from early_pole_decay_mode_clock_reference import G, KAPPA, intervals


def bare_closed(C, f, E, m, beta):
    M,b=.05,.003
    ep=np.hypot(C.M.edge_sets[1],M)
    eb=np.hypot(C.M.edge_sets[2],b)
    z,w=roots_legendre(16)
    vals=[]
    for q,es,N,V in zip(C.q,E,C.ns,C.sw,strict=True):
        _,l,h,*_=intervals(M,m,b,q)
        edges=np.unique(np.r_[l,h,ep[(ep>l)&(ep<h)],
                             eb[(eb+es>l)&(eb+es<h)]+es])
        half=np.diff(edges)[:,None]/2
        x=((edges[1:]+edges[:-1])[:,None]/2+half*z).ravel()
        ww=(half*w).ravel()
        def shell(en,ee,j):
            i=np.searchsorted(ee,en)-1
            k=C.M.r[j].start+np.clip(i,0,len(ee)-2)
            u=(en-C.M.shape[k,0])/C.M.shape[k,1]
            return np.where((i>=0)&(i<len(ee)-1),expit(f[k,0]+f[k,1]*u),0.)
        F, B=shell(x,ep,1),shell(x-es,eb,2)
        vals.append(beta*(N*(ww@(F-B))+V*(ww@(F*(1-B))))/(q*es))
    return float(sum(vals))


def evaluate(profile,subdivision):
    grid=profile['grid']
    with gzip.open(GRID/f'early-pole-spectrum-{grid}-m06-2026-11-26.json.gz','rt') as f:
        a=json.load(f)
    archive=np.asarray(a['endpoint_parent_momentum_edges'])
    f=np.asarray(a['endpoint_parent_occupations'])
    e,p,w,_=refine(archive,f,profile['parent'],a['x_end'])
    C=Cold(grid,1,0,8,True,aperture_subdivision=subdivision)
    Y0=C.M.init.copy()
    selection=C.M.r[1]
    Y1=Y0.copy()
    ids=[]
    for i in range(selection.start,selection.stop):
        K=C.M.cells[i]
        if not any(K.lo<t['qEdges'][-1] and K.hi>t['qEdges'][0] for t in profile['parent']):
            continue
        old=fit_archive_cell(K,e,p)
        np.testing.assert_allclose(old,Y0[i],rtol=3e-10,atol=2e-16)
        Y1[i]=fit_archive_cell(K,e,w)
        ids.append(i)
    E,_,_,_,m=C.scalar(C.ns)
    beta=G*KAPPA*.05*((1.06*.05)**2-m*m)
    d_edges=np.asarray(a['endpoint_fermion_momentum_edges'])
    d_f=np.asarray(a['endpoint_fermion_occupations'])
    def bare_raw(alpha):
        v=[]
        for q,N,V in zip(C.q,C.ns,C.sw,strict=True):
            es,I,_=mode(.05,m,.003,q,e,p,w,d_edges,d_f,alpha)
            v.append(beta*(N*(I[0]-I[1])+V*I[2])/(q*es))
        return float(sum(v))
    def audit(Y,alpha):
        e,f,grad,rho,S=C.M.weights(Y,True)
        k=np.r_[Y.ravel(),C.ns,0.,0.]
        rate=Cold.rhs(C,k,False)
        d=rate[C.j0:C.j0+C.n]
        return dict(hard_energy=rho,hard_entropy=S,
                    cold_rate=float(d.sum()),cold_rate_by_q=d.tolist(),entropy_rate=rate[-2],
                    bare_projected_cold_rate=bare_closed(C,f,E,m,beta),
                    bare_raw_cold_rate=bare_raw(alpha),
                    forward_events_energy_defect=rate[-1])
    branches=[dict(witness_moment_mix=alpha,**audit((1-alpha)*Y0+alpha*Y1,alpha)) for alpha in (0.,.5,.85,.97)]
    assert abs(np.sum((Y1-Y0)[selection,0]))<2e-15
    # Compare integrated bare energy in the two moment states by exact
    # affine representation within each live cell, before dressing.
    bare=np.asarray([[C.M.cells[i].ec,C.M.cells[i].de] for i in ids])
    diff=float(np.sum(bare*(Y1-Y0)[ids]))
    assert abs(diff)<2e-15
    return dict(grid=grid,subdivision=subdivision,changed_live_parent_cells=len(ids),
                parent_count_defect=float(np.sum((Y1-Y0)[selection,0])),bare_E_defect=diff,
                q=C.q.tolist(),branches=branches)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=('fine','third'))
    p.add_argument('--subdivisions',type=int,nargs='+',default=[1,2,4,8])
    a=p.parse_args()
    profile=next(p for p in json.loads(PROFILES.read_text()) if p['grid']==a.grid)
    print(json.dumps([evaluate(profile,i) for i in a.subdivisions],indent=2,allow_nan=False))
