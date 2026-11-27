"""Local exact linear ramp in phase-space volume, without artificial subcell edges."""
import argparse
import gzip
import json

import numpy as np
from scipy.special import roots_legendre

from early_pole_archive_scalar_below_pole_finite_q_reference import integrate
from early_pole_archive_scalar_pole_core_reference import energy_integral,residual,fixture
from early_pole_native_aperture_witness_reference import GRID,occup


def exact_ramp(a,g,frac):
    p=np.asarray(a['endpoint_fermion_momentum_edges'])
    f=np.asarray(a['endpoint_fermion_occupations'])
    edge=1.4222222222222223 if g=='fine' else 1.44
    i=np.argmin(abs(p-edge))
    assert abs(p[i]-edge)<1e-10
    v=p**3
    dx=frac*min(v[i]-v[i-1],v[i+1]-v[i])
    left,right=v[i]-dx,v[i]+dx
    ep=np.r_[v[:i],(v[i-1]+left)/2,left,right,(right+v[i+1])/2,v[i+1:]]
    out=np.r_[f[:i],f[i-1],(f[i-1]+f[i])/2,f[i],f[i:]]
    dN=np.diff(ep)
    dE=np.diff(energy_integral(ep))
    x,w=roots_legendre(16)
    correction=dx*((f[i]-f[i-1])*.5*x*np.hypot(np.cbrt(v[i]+dx*x),.003))@w
    target=np.array([np.diff(v[i-1:i+2])@f[i-1:i+1],
                     np.diff(energy_integral(v[i-1:i+2]))@f[i-1:i+1]])
    slots=slice(i-1,i+4)
    ids=[i-1,i+3]
    delta=np.linalg.solve(np.array([dN[ids],dE[ids]]),target-np.array([dN[slots]@out[slots],dE[slots]@out[slots]+correction]))
    out[ids]+=delta
    np.testing.assert_allclose([dN[slots]@out[slots],dE[slots]@out[slots]+correction],target,rtol=5e-14)
    assert np.all((out>=0)&(out<=1))
    pe=np.cbrt(ep)

    def fermi(p):
        v=p**3
        r=(v-left)/(right-left)
        return np.where((left<=v)&(v<right), (1-r)*f[i-1]+r*f[i],occup(p,pe,out))
    trial=dict(a,endpoint_fermion_momentum_edges=pe.tolist(),endpoint_fermion_occupations=out.tolist())
    return trial,fermi,dict(ramp_momentum_interval=np.cbrt([left,right]).tolist(),
                            occupation_adjustments_on_shoulders=delta.tolist(),moment_weights=target.tolist())


def run(g,frac):
    with gzip.open(GRID/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz','rt') as fh:
        a=json.load(fh)
    a,fermi,info=exact_ramp(a,g,frac)
    q=4.1e-6
    M,m=.05,.003
    ms=a['checkpoints'][-1]['scalar_pole_over_T']
    ps=np.array(a['endpoint_fermion_momentum_edges'])
    es=np.hypot(ps,m)
    radius=.002
    sides=np.abs((m*m+ms*ms+2*(es*np.hypot(q,ms)+np.array([-1,1])[:,None]*ps*q))/(M*M)-1)
    breaks=np.unique(sides[sides<radius])
    out=[]
    for order in (18,28):
        def fetch(x):
            r=integrate(a,q,order,.9,probe_x=x,resolve_terms=True,fermi_shape=fermi)
            return np.concatenate([np.array(r['term_gross_densities_per_reference_H']),
                                   np.array(r['term_net_densities_per_reference_H'])[:,:,None]],axis=2)
        memo={}
        def get(x):
            if x not in memo:
                memo[x]=fetch(x)
            return memo[x]
        F=get(0.)
        even=[]
        for x in (2e-4,1e-4,2e-5,1e-5,2e-6,1e-6,2e-7,2e-8):
            v=(get(x)[0]-F[0])+(get(-x)[0]-F[0])
            even.append(dict(x=x,even_over_x=(v/x).tolist(),even_over_x2=(v/x**2).tolist()))
        for r in (.0002,radius):
            for n in (8,14):
                for ratio in (0.,.1):
                    v=residual(r,ratio,n,F,get,breaks)
                    out.append(dict(order=order,window=r,core_order=n,ratio=ratio,terms=v.tolist(),total=v.sum(0).tolist()))
        info[f'order{order}_parity']=even
    return dict(grid=g,ramp_fraction=frac,**info,topology_thresholds=breaks.tolist(),cases=out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=('fine','third'));p.add_argument('--fraction',type=float,default=.03)
    a=p.parse_args()
    if not 0.<a.fraction<.5:
        p.error("The ramp fraction must be strictly between 0 and .5.")
    fixture()
    print(json.dumps(run(a.grid,a.fraction),indent=2,allow_nan=False))
