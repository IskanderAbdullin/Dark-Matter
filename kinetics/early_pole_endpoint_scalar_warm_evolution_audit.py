"""Compare matched cold and warm-only refined scalar histories at x=.05."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_scalar_spectral_moment_budget import LABELS,integrate
from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_scalar_shape_transport import CUT,cdf,pair
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart
from early_pole_endpoint_stencil_remap_pilot import SPECS,X,load


def unpack(c,warm,x=.05):
    if (c['grid'],c['charge_div'],c['scalar_div'],c['cut'],c['reference_x'],c['x_end'],c['steps'])!=(
            'fine',2,4,1e-5,X,x,1024 if x==.01 else 2048):
        raise ArithmeticError('Different cold, charge or time setting')
    if (not warm and 'warm_floor' in c) or (warm and (c.get('warm_floor')!=.01 or c.get('warm_div',0)<2)):
        raise ArithmeticError('Expected an unrefined warm chart followed by a warm-only refinement')
    opt={k:c[k] for k in ('warm_floor','warm_div')} if warm else {}
    d=ScalarSubshellChart(load('fine')[0],SPECS['fine'],2,4,1e-5,**opt)
    a=np.asarray(c['endpoint_occupation'])
    if c['counts']!=[d.ns,d.nf,d.np] or a.shape!=(d.n,) or not np.all(np.isfinite(a)) or min(a)<0 or max(a[d.ns:])>1:
        raise ArithmeticError('Different state from declared endpoint chart')
    return d,a


def warm_cdf(d,y,lo=.01):
    native=d.ds.ed
    if not native[0]<lo<native[-1]:
        raise ArithmeticError('Warm cut lies outside the retained scalar archive')
    e=np.r_[lo,native[native>lo]]
    l,r=e[:-1],e[1:]
    j=np.searchsorted(native,(l+r)/2)-1
    n=y[j]*(r-l)*(r*r+r*l+l*l)
    M=fsum(n)
    if not M>0:
        raise ArithmeticError('No warm population to normalize')
    u=n/M
    mean=.75*fsum(y[j]*(r-l)*(r**3+r*r*l+r*l*l+l**3))/M
    return e,np.r_[0.,np.cumsum(u)],np.r_[np.cumsum(u[::-1])[::-1],0.],mean


def run(v,w,x=.05):
    (d,a),(e,b)=unpack(v,False,x),unpack(w,True,x)
    if any(not np.array_equal(getattr(d,k),getattr(e,k)) for k in ('qp','qd','wp','wd')) or (
            not np.array_equal(d.ds.ed[d.ds.ed<=.01],e.ds.ed[e.ds.ed<=.01])):
        raise ArithmeticError('Refinement changed fixed charge nodes or the cold grid')
    m=d.state(x,a)[0]
    pairs=((d,a),(e,b))
    u,z=[integrate(h,s,m,48) for h,s in pairs]
    for (h,s),ref in zip(pairs,(u,z),strict=True):
        check=integrate(h,s,m,24)
        if np.max(abs(check-ref)/np.maximum(abs(ref),1e-105))>2e-12:
            raise ArithmeticError('An actual warm moment is unresolved')
    exact=np.stack((u.sum(axis=0),z.sum(axis=0)))
    def nodal(h,s):
        q=h.ds.p
        E=np.hypot(q,m)
        N=h.ds.w*s[:h.ns]
        X=np.stack((q**0,q,q*q,q*q/(E+m),q*q/E),axis=-1)
        return np.array([fsum(N*X[:,j]) for j in range(len(LABELS))])
    node=np.array([nodal(h,s) for h,s in pairs])
    if np.max(abs(node[:,0]/exact[:,0]-1))>2e-14:
        raise ArithmeticError('Exact and nodal count did not agree')
    return dict(x=x,warm_div=w['warm_div'],warm_floor=w['warm_floor'],fixed_mass=float(m),
                scalar_counts=[d.ns,e.ns],columns=LABELS,
                cold_conditional_transport=pair(cdf(d,a,CUT),cdf(e,b,CUT)),
                warm_conditional_transport=pair(warm_cdf(d,a),warm_cdf(e,b)),
                total_relative_changes=(exact[1]/exact[0]-1).tolist(),
                nodal_total_relative_changes=(node[1]/node[0]-1).tolist(),
                base_and_refined_nodal_over_exact_minus_one=(node/exact-1).tolist(),
                relative_changes_at_q_at_least_point01=(z[-2:].sum(axis=0)/u[-2:].sum(axis=0)-1).tolist(),
                scope='Warm-only evolution sensitivity on a fixed cold and charged chart. No physical convergence claim.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('s4',type=Path)
    p.add_argument('s4_warm',type=Path)
    p.add_argument('--x',type=float,choices=(.01,.02,.04,.05),default=.05)
    a=p.parse_args()
    print(json.dumps(run(read(a.s4),read(a.s4_warm),a.x),indent=2,allow_nan=False))
