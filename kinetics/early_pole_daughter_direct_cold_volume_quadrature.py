"""Integrate the raw cold seed density in scalar volume on held native packets."""
import argparse
import gzip
import json
from math import fsum
from types import SimpleNamespace

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from scipy.special import expit, roots_legendre

from early_pole_decay_mode_clock_reference import G,KAPPA
from early_pole_daughter_local_fraction_curvature_attribution import initial
from early_pole_frozen_dynamic_vlasov_reference import load
from early_pole_live_pauli_gauge_moment_decay_reference import Tool


def angular_edge_modes(M,e,mass,root_scan):
    lo,hi=1e-6,1e-5
    q=np.linspace(lo,hi,root_scan)
    knots=[]
    for i in range(M.r[1].start,M.r[1].stop):
        A=M.cells[i]
        ea=lambda p:e[i,0]+e[i,1]*A.basis(p)
        for k in range(M.r[2].start,M.r[2].stop):
            B=M.cells[k]
            if (B.lo<1.4138888888888888-2e-14 or B.hi>1.4666666666666666+2e-14
                    or max(A.lo-B.hi,B.lo-A.hi)>hi):
                continue
            eb=lambda p:e[k,0]+e[k,1]*B.basis(p)
            for fixed,other,fa,fb,direction in ((A,B,ea,eb,1.),(B,A,eb,ea,-1.)):
                for endpoint in (fixed.lo,fixed.hi):
                    for sign in (-1.,1.):
                        test=lambda x:direction*(fa(endpoint)-fb(endpoint+sign*x))-np.hypot(mass,x)
                        v=test(q)
                        for j in np.flatnonzero(v[:-1]*v[1:]<0.):
                            root=brentq(test,q[j],q[j+1],xtol=1e-17)
                            p=endpoint+sign*root
                            if other.lo-1e-13<=p<=other.hi+1e-13:
                                knots.append(root)
    # Distinct geometric cuts that are within numerical boundary accuracy are merged.
    value=np.sort(np.asarray(knots))
    if not len(value):
        return value
    return value[np.r_[True,np.diff(value)>1e-12]]


def direct_cold(archive,order,sub,split_angular_edges=False,root_scan=129,shared_angular_edges=None):
    C,seed,(Y,_,_),_,depth=initial(archive)
    M=C.M
    _,_,_,_,mass=C.scalar(seed)
    e,F=M.weights(Y,True)[:2]
    with gzip.open(archive,'rt') as f:
        grid=json.load(f)['branches']['steps']['grid']
    d=load(f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz')
    p=np.asarray(d['endpoint_fermion_momentum_edges'])
    bf=np.asarray(d['endpoint_fermion_occupations'])
    p3=p**3/3.
    u=PchipInterpolator(p3,np.r_[0.,np.cumsum(bf*np.diff(p3))],extrapolate=False).derivative()
    cd=np.asarray(d['scalar_momentum_edges'])**3
    nn=np.maximum(d['endpoint_scalar_signed_shell_number_over_T3'],0.)
    if nn.shape!=(len(cd)-1,):
        raise ValueError('Incompatible scalar archive support.')
    if shared_angular_edges is not None:
        knots=np.asarray(shared_angular_edges)
        if (knots.ndim!=1 or not np.all(np.isfinite(knots)&(knots>1e-6)&(knots<1e-5))
                or not np.all(np.diff(knots)>0.)):
            raise ValueError('Invalid shared angular cuts.')
    else:
        knots=angular_edge_modes(M,e,mass,root_scan) if split_angular_edges else np.array([])
    v=np.unique(np.r_[np.geomspace(1e-6,1e-5,9)**3,cd[(cd>1e-18)&(cd<1e-15)],knots**3])
    # Scalar number weight is d(q^3)/(6*pi^2).
    x=.05*np.exp(getattr(M,"eta",0.))
    pref=G*KAPPA*x*((1.06*x)**2-mass*mass)/(6*np.pi**2)
    node,weight=roots_legendre(order)
    initial_groups=np.zeros(8)
    integrated=[[] for _ in range(8)]
    q8=np.geomspace(1e-6,1e-5,9)**3
    term_count=0
    for a,b in zip(v[:-1],v[1:],strict=True):
        m=(a+b)/2
        j=np.searchsorted(cd,m)-1
        shell=np.searchsorted(q8,m)-1
        if not (0<=j<len(nn) and 0<=shell<8):
            raise ValueError('Volume outside two matched archives.')
        density=nn[j]/(cd[j+1]-cd[j])
        scalar_f=6*np.pi**2*density
        initial_groups[shell]+=(b-a)*density
        subgrid=np.linspace(a,b,sub+1)
        for lo,hi in zip(subgrid[:-1],subgrid[1:],strict=True):
            dv=(hi-lo)/2
            for V,w in zip((lo+hi)/2+dv*node,dv*weight,strict=True):
                q=np.cbrt(V)
                deg=np.hypot(q,mass)
                model=SimpleNamespace(n=M.decay.n,cells=M.decay.cells,q=q,nod=M.z,w=M.w)
                for i,k,xi,xf,ww in Tool.events(model,dict(coeff=e[M.start:],es=deg)):
                    i+=M.start
                    k+=M.start
                    z=M.cells[k]
                    if z.lo<1.4138888888888888-2e-14 or z.hi>1.4666666666666666+2e-14:
                        continue
                    pp=z.pfrombare(z.ec+z.de*xf)
                    ids=np.searchsorted(p,pp)-1
                    B=np.maximum(u(pp**3/3.),0.)-bf[ids]
                    if np.any(~np.isfinite(B)) or np.any((ids<0)|(ids>=len(bf))):
                        raise ValueError('Out-of-support daughter quadrature.')
                    Pa=expit(F[i,0]+F[i,1]*xi)
                    term=-pref*w/(q*deg)*float(ww@((Pa+scalar_f)*B))
                    integrated[shell].append(term)
                    term_count+=1
    np.testing.assert_allclose(initial_groups,seed,rtol=0.,atol=5e-18)
    sums=[fsum(g) for g in integrated]
    return dict(**(dict(shared_scalar_volume_edges=v.tolist(),scalar_evaluations=(len(v)-1)*order*sub)
                   if shared_angular_edges is not None else {}),
                hard_depth=depth,scalar_order=order,subpanels_per_base_cut=sub,
                angular_edge_modes=knots.tolist(),angular_root_scan=root_scan if split_angular_edges else None,
                initial_pole=float(mass),grouped_seed_error=float(max(abs(initial_groups-seed))),
                direct_scalar_volume_rate=fsum(sums),grouped_rate=sums,packet_cold_terms=term_count)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('archive',nargs='+')
    p.add_argument('--order',type=int,default=8)
    p.add_argument('--subpanels',type=int,default=1)
    p.add_argument('--split-angular-edges',action='store_true')
    p.add_argument('--root-scan',type=int,default=129)
    v=p.parse_args()
    if v.order<2 or v.subpanels<1 or v.root_scan<33:
        p.error('Require positive scalar order and panel count.')
    print(json.dumps([direct_cold(f,v.order,v.subpanels,v.split_angular_edges,v.root_scan) for f in v.archive],indent=2,allow_nan=False))
